# this script first generate Amphibian and Reptile AOH 
# first extract land use info to spp shapefile 
# then extract DEM info to spp shapefile 
# cut to areas overlapping with forest wetland cropland 
# further contrain to areas with elevation 0 - 1700m 
# then calculate area in resolution of 1000 m

library(sf)
library(rgee)
library(terra)
library(tidyverse)

ee_Initialize()
# ee_Authenticate() 
ee_check()

# read data 
species.list <- list.files("./data/trade_amphi", pattern = ".shp")


# get land use data
landuse <- ee$ImageCollection('ESA/WorldCover/v200')$first()
# urban accessibility data
urban_accessibility <- ee$Image("Oxford/MAP/accessibility_to_cities_2015_v1_0")
# DEM data
srtm <- ee$Image("CGIAR/SRTM90_V4")

# some input polygon is too big for extraction. to avoid memory issues, split into 5x5 degree grids 
# Define the grid size (e.g., 1 degree grid or 10km grid in meters)
grid_size <- 5  # Adjust the grid size to control chunk size (degrees for EPSG:4326)

# # Check the split polygons
# plot(split_polygons, col = "lightblue", border = "white")

#AOH.list <- tibble()
for (i in (337:length(species.list))) {
  # read species shapefile 
  spp <-  st_read(paste0("./data/trade_amphi/",species.list[i]))
  
  # Create a grid over the bounding box of the polygon
  grid <- st_make_grid(spp, cellsize = grid_size, what = "polygons")
  # Intersect the grid with the original polygon to get smaller parts
  split_spp <- st_intersection(grid, spp)
  
  AOH.spp.i = tibble()
  for (ii in (1:length(split_spp))) {
    print(paste0("processing spp ", i, ", part ", ii, " out of ", length(split_spp), " parts."))
    # turn to gee file 
    spp_ee <- sf_as_ee(split_spp[ii])
    # convert to raster to have every pixel with an initial value of 0
    spp_r <- ee$Image$constant(0)$clip(spp_ee)$reproject(
      crs = "EPSG:4326",
      scale = 1000
    )
    
    # ---- land use ---- #
    # Create a mask for land cover values 90, 40, or 10
    landcover_mask <- landuse$eq(90)$Or(landuse$eq(40))$Or(landuse$eq(10))
    # Update spp.r where the mask is true
    spp_r <- spp_r$where(landcover_mask, 1)
    
    # ---- DEM --- #
    # Create a mask for elevation between 0 and 1700
    elevation_mask <- srtm$gte(0)$And(srtm$lte(1700))
    # Update spp.r: keep only pixels where elevation mask is true
    spp_r <- spp_r$updateMask(elevation_mask)
    
    # Calculate the number of pixels equal to 1
    pixel_count <- spp_r$eq(1)$reduceRegion(
      reducer = ee$Reducer$sum(),
      geometry = spp_ee,
      scale = 1000,
      maxPixels = 1e13
    )$get("constant")  # "constant" refers to the value in the raster
    
    # Fetch the result
    pixel_count_result <- pixel_count$getInfo()
    
    ## ------ urban accessibility ------ ##
    masked_accessibility <- urban_accessibility$updateMask(spp_r$eq(1))
    # Compute the mean and standard deviation of urban accessibility
    stats <- masked_accessibility$reduceRegion(
      reducer = ee$Reducer$mean()$combine(
        reducer2 = ee$Reducer$stdDev(),
        sharedInputs = TRUE
      ),
      geometry = spp_ee, 
      scale = 1000,                    # Match the raster's resolution
      maxPixels = 1e13                 # Handle large datasets
    )
    
    # Fetch the results
    mean_access <- stats$get("accessibility_mean")$getInfo()          # Mean of urban accessibility
    stddev_access <- stats$get("accessibility_stdDev")$getInfo()      # Standard deviation of urban accessibility

    if (is.null(mean_access)) {
      AOH.spp.ii = tibble(i = i,
                       ii = ii,
                       id_no = spp$id_no, sci_name = spp$sci_name, 
                       AOH_count = round(pixel_count_result,2),
                       mean_access = NA,
                       stddev_access = NA)
    } else {
      AOH.spp.ii = tibble(i = i,
                       ii = ii,
                       id_no = spp$id_no, 
                       sci_name = spp$sci_name, 
                       AOH_count = round(pixel_count_result,2),
                       mean_access = mean_access,
                       stddev_access = stddev_access)
    }
     
    AOH.spp.i = rbind(AOH.spp.i, AOH.spp.ii) 
    }

  AOH.list = rbind(AOH.list, AOH.spp.i)
}

# AOH.list.first <- read_rds("./results/AOH_up_tp_45.rds") %>% distinct() %>% mutate(i = row_number(.), ii = 0)
# AOH.list.full = rbind(AOH.list.first, AOH.list %>% distinct())



# now calculate overall mean and stddev access for each species 
stats <- AOH.list %>% 
  group_by(sci_name) %>% 
  summarise(
    # Step 1: Weighted mean
    spp_access_mean = sum(AOH_count * mean_access) / sum(AOH_count),
    
    # Step 2: Within-group variance
    within_var = sum((AOH_count - 1) * stddev_access^2) / (sum(AOH_count) - n()),
    
    # Step 3: Between-group variance
    between_var = sum(AOH_count * (mean_access - (sum(AOH_count * mean_access) / sum(AOH_count)))^2) / sum(AOH_count),
    
    # Step 4: Population variance
    spp_access_var = within_var + between_var,
    
    # Step 5: Population standard deviation
    spp_access_sd = sqrt(spp_access_var),
    
    # Total sample size
    total_AOH_count = sum(AOH_count),
    .groups = "drop"
  )

write_rds(AOH.list, "./results/AOH_amphi_allspp.rds")
write_csv(stats %>% dplyr::select(sci_name, total_AOH_count, spp_access_mean, spp_access_var), "./results/trade_amphi_results.csv")
