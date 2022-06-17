#### Title: "Converting GAO imagery to larval cloud"
#### author: "Rachel Carlson"
#### date: "6 February, 2022"
#### summary: The purpose of this script is to convert GAO data for my tile region to a point cloud for larval modeling in Ocean Parcels.


library(raster)
library(sf)
library(tidyverse)

## Import GAO data

# All coral genera
gao <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/ASU_GAO_Hawaii_SW_Live_v2.tif")
gao_W <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/GAO input/West_Hawaii_merged_CC.tif")
ncell(gao) # There are 1,075,206,200 cells in the SW Hawaii raster and 2,247,402,650 cells in the NW Hawaii raster (resolution 2 m)
# Porites for SW Hawaii
# Pocillopora for SW Hawaii
# Montipora for SW Hawaii

## Resample GAO data from 2 m to 40 m, 100 m, and 200 m
gao_40m <- aggregate(gao, fact = 20)
gao_100m <- aggregate(gao, fact = 50)
gao_200m <- aggregate(gao, fact = 100)

RastertoSeed <- function(raster) {
  # Function to convert raster to an array of xy coordinates (in m) representing the center of each pixel and pixel value.
  points <- rastertoPoints(raster) ## Convert raster pixels to point cloud
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
  return(sf_larv)
}


## Export to shapefile
st_write(sf_40m, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud.shp")

## Coordinates are in meters, so need to redefine as decimal degrees for Parcels
## From importing the point shapefile above into QGIS, we know the original coordinate system = 32605
st_crs(sf_40m) = 32605 # We need to define this before reprojecting (otherwise won't know how to reproject)

## Then reproject to a coordinate system with decimal degrees (4326)
cloud_40 <- st_transform(sf_40m, 4326)
View(cloud_40) # It appears to have worked. View in QGIS to check (it does).

## Finally, since Parcels needs lat and lon in separate fields, define lat and long columns from geometry
cloud_40 <- cloud_40 %>% mutate(lat = unlist(map(cloud_40$geometry,2)),
         lon = unlist(map(cloud_40$geometry,1)))

## This reduces lat/lon to 4 significant figures, but upon viewing in QGIS, this retains the grid structure and is accurate to 40m

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

## Read in one layer of current data and convert to one big polygon (HYCOM ocean boundary)
foo <- raster("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/oceanographic/HYCOM/Kona_uvwts_300.cdf")
pp <- rasterToPolygons(foo, dissolve = TRUE)
pfoo <- pp %>% st_as_sf()
pfoo2 <- st_union(pfoo)
plot(pfoo2)

## Intersect points and HYCOM boundary polygon
sf2 <- st_intersects(sf, pfoo2) # This gives you a list of 15600 points where row is 1 if intersects, (empty) if doesn't intersect
b <- sapply(sf2,function(x){length(x)==0}) # This gives a logical vector of length 15600 with TRUE/FALSE intersects
# c <- !b # Use ! because we want the inverse of this vector to find those points that DO NOT intersect
sf3 <- sf[b,] # Subset our points by only those that are outside the HYCOM ocean boundary

## Snap those points to the nearest HYCOM boundary line
# This didn't work -- use ArcGIS `Near` tool following directions here https://support.esri.com/en/technical-article/000021426 to find XNear/YNear (new x and y for points outside of HYCOM)
snapped <- st_read("/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/gao_pointcloud_180_snapped.shp") %>% 
  as_data_frame() %>% 
  select(-c(geometry,NEAR_FID)) %>% # Default geometry is original lat/lon, need to reset
  st_as_sf(coords = c("NEAR_X", "NEAR_Y")) %>% # Convert to CSV and back to sf with new, nearest HYCOM boundary lat/lon
  sf::st_set_crs(4326)
# This only represents those points outside the boundary originally

## Join these points to larger dataset
sf4 <- sf[!b,]
GAO_cloud <- bind_rows(sf4,snapped) %>% select(-c(lat,lon))
sum(!is.na(GAO_cloud$NEAR_DIST)) # Check to make sure this is the # of snapped points (e.g., 7110)
st_write(GAO_cloud, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud/gao_plus_snapped.shp")


