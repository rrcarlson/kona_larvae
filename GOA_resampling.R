#### Title: "Converting GAO imagery to larval cloud"
#### author: "Rachel Carlson"
#### date: "6 February, 2022"
#### summary: The purpose of this script is to convert GAO data for my tile region to a point cloud for larval modeling in Ocean Parcels.


library(raster)
library(sf)
library(tidyverse)

## Import GAO data

# All coral genera. Added aggregate(fact = 20) to resample GAO to 40 m
gao <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/ASU_GAO_Hawaii_SW_Live_v2.tif") %>% aggregate(fact = 20)
gao_W <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/West_Hawaii_merged_CC.tif") %>% aggregate(fact = 20)
ncell(gao) # There are 1,075,206,200 cells in the SW Hawaii raster and 2,247,402,650 cells in the NW Hawaii raster (resolution 2 m)
# Porites for SW Hawaii
Por <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/ASU_GAO_SfMCover_v1_All_Porites.tif") %>% aggregate(fact = 20)
# Pocillopora for SW Hawaii
Poc <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/ASU_GAO_SfMCover_v1_All_Pocillopora.tif") %>% aggregate(fact = 20)
# Montipora for SW Hawaii
Mon <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/ASU_GAO_SfMCover_v1_All_Montipora.tif") %>% aggregate(fact = 20)

# gao_100m <- aggregate(gao, fact = 50)
# gao_200m <- aggregate(gao, fact = 100)

## Create sf object for HYCOM boundary (where current vectors are/not NA)
# Read in one layer of current data and convert to one big polygon
HYCOM_bound <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/oceanographic/HYCOM/Kona_uvwts_300.cdf")
pp <- rasterToPolygons(HYCOM_bound, dissolve = TRUE)
pf <- pp %>% st_as_sf() %>% st_union()
plot(pf) # Solid, single polygon showing HYCOM boundary

################ Functions

### Function to convert GAO raster to larval point cloud (sf_larv, an sf object)
RastertoSeed <- function(raster) {
  # Function to convert raster to an array of xy coordinates (in m) representing the center of each pixel and pixel value.
  points <- raster::rasterToPoints(raster) ## Convert raster pixels to point cloud
  sf_larv <- points %>% as.data.frame %>% sf::st_as_sf(coords = c(1,2)) # Function to convert array to sf object (point shapefile)
  colnames(sf_larv) <- c("live_coral", "geometry")
  # Create attribute for number of coral larvae each point represents
          ## The dataframe above has one attribute per point: % live coral cover, ranging from 0 - 64.45%.
          ## Larvae at each point is proportional to coral cover. Since there are no fractional larvae,
          ## I will first translate coral cover (i.e., larval amount) to the nearest integer.
          ## Then, for computational efficiency, I will use each larvae to represent ten larvae, so will find the
          ## nearest integer to 1/10 coral cover.
  sf_larv$larvae <- round(sf_larv$live_coral) # Translate coral cover to nearest integer.
  sf_larv$larv_red <- round(sf_larv$larvae/10) # Same as above, but if all larvae represent 10 larvae (faster computation)
  ## Coordinates are in meters, so need to redefine as decimal degrees for Parcels
  ## From importing the point shapefile above into QGIS, we know the original coordinate system = 32605
  st_crs(sf_larv) <- 32605 # We need to define this before reprojecting (otherwise won't know how to reproject)
  sf_larv <- st_transform(sf_larv, 4326) # Reproject to a coordinate system with decimal degrees (4326)
  sf_larv <- sf_larv %>% mutate(lat = unlist(map(sf_larv$geometry,2)), # Since Parcels needs lat and lon in separate fields, define lat and lon columns from geometry
                                lon = unlist(map(sf_larv$geometry,1))) # This reduces lat/lon to 4 significant figures but it's OK; in QGIS, we can see that the correct point cloud structure is retained and accurate to 40m
  sf_larv$newlon <- sf_larv$lon + 360 # We'll use -180/180 version here, but this gives us the option of using 0/360
  return(sf_larv) 
}

foo <- RastertoSeed(gao_40m)


### Function to subset larvae points that fall in NA area of HYCOM data
# hycom_bound = a polygon representing the boundary of current data
# sf_larv = an sf object representing a larval cloud (output of RastertoSeed)
# snapped = fp for a shapefile processed from sf_larv in ArcGIS using the 'Near' tool. These are points from output of RastertoSeed, snapped to the nearest HYCOM boundary line (follow https://support.esri.com/en/technical-article/000021426 to find XNear/YNear, new x and y for points outside of HYCOM)
PointsOutside <- function(hycom_bound, sf_larv, snapped_fp) {
  # Intersect points and HYCOM boundary polygon
  inter_points <- st_intersects(sf_larv, hycom_bound) # This gives you a list of 15600 points (for all coral) where row is 1 if intersects, (empty) if doesn't intersect
  b <- sapply(inter_points,function(x){length(x)==0}) # This gives a logical vector of length 15600 (or other sf_larv length) with TRUE/FALSE intersects
  sf_in <- sf_larv[!b,] # !b is points inside HYCOM boundary. !b SHOULD be the inverse of the b vector, i.e., those points that DO NOT intersect-but weirdly, it's the opposite
  snapped <- st_read(snapped_fp) %>% # Snapped will end up only being points outside the HYCOM boundary
    as_data_frame() %>% 
    select(-c(geometry,NEAR_FID)) %>% # Default geometry is original lat/lon, need to reset
    st_as_sf(coords = c("NEAR_X", "NEAR_Y")) %>% # Convert to CSV and back to sf with new, nearest HYCOM boundary lat/lon
    sf::st_set_crs(4326)
  final_cloud <- bind_rows(sf_in,snapped) %>% select(-c(lat,lon))
  return(final_cloud)
  }


# QA/QC
sum(!is.na(final_cloud$NEAR_DIST)) # Check to make sure this is the # of snapped points (e.g., 7110)
nrow(snapped) # Should be identical
final_cloud_csv <- as.data.frame(cloud_40) %>% select(-geometry) # Finally, convert to a csv without geometry column for use in Parcels

st_write(final_cloud, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/gao_plus_snapped.shp")
write.csv(final_cloud_csv, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud.csv", row.names = FALSE) # row.names = FALSE removes unlabeled index column in CSV






























## Export to shapefile
st_write(sf_40m, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud.shp")

## Coordinates are in meters, so need to redefine as decimal degrees for Parcels
## From importing the point shapefile above into QGIS, we know the original coordinate system = 32605
st_crs(sf_40m) = 32605 # We need to define this before reprojecting (otherwise won't know how to reproject)

## Then reproject to a coordinate system with decimal degrees (4326)
cloud_40 <- st_transform(sf_40m, 4326)
View(cloud_40) # It appears to have worked. View in QGIS to check (it does).

## Since Parcels needs lat and lon in separate fields, define lat and lon columns from geometry
cloud_40 <- cloud_40 %>% mutate(lat = unlist(map(cloud_40$geometry,2)),
                                lon = unlist(map(cloud_40$geometry,1)))

## This reduces lat/lon to 4 significant figures but it's OK; in QGIS, we can see that the correct point cloud structure is retained and accurate to 40m

## Finally, convert to a csv without geometry column for use in Parcels
df <- as.data.frame(cloud_40)[, c(1, 3:6)]
write.csv(df, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud.csv", row.names = FALSE) # row.names = FALSE removes unlabeled index column in CSV

#################### Snap point cloud to non-NA current data

## Read in point cloud
df <- read.csv("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud.csv")
df$newlon <- df$lon + 360 # We'll use -180/180 version here, but this gives us the option of using 0/360

## Convert to point sf object
sf <- df %>% st_as_sf(coords = c("lon", "lat")) %>%
  sf::st_set_crs(4326)