#### Title: "Converting GAO imagery to larval cloud"
#### author: "Rachel Carlson"
#### date: "6 February, 2022"
#### summary: The purpose of this script is to convert GAO data for my tile region to a point cloud for larval modeling in Ocean Parcels.


library(raster)
library(sf)
library(tidyverse)

## Import GAO data
gao <- raster("/Users/rachelcarlson/Documents/Research/RS_data/West_Hawaii_CC/West_Hawaii_Island/LiveCoral/ASU_GAO_Hawaii_SW_Live_v2.tif")
ncell(gao) # There are 1,075,206,200 cells in this raster (resolution 2 m)

## Resample GAO data from 2 m to 40 m, 100 m, and 200 m
gao_40m <- aggregate(gao, fact = 20)
gao_100m <- aggregate(gao, fact = 50)
gao_200m <- aggregate(gao, fact = 100)
writeRaster(gao_40m, "/Users/rachelcarlson/Documents/Research/RS_data/West_Hawaii_CC/West_Hawaii_Island/LiveCoral/40m_ASU_GAO_Hawaii_SW_Live_v2.tif")

## Convert raster pixels to point cloud
points_40m <- rasterToPoints(gao_40m) # Function to convert raster to an array of xy coordinates (in m) representing the center of each pixel and pixel value.
sf_40m <- points_40m %>% as.data.frame %>% sf::st_as_sf(coords = c(1,2)) # Function to convert array to sf object (point shapefile)
colnames(sf_40m) <- c("live_coral", "geometry")

## Create attribute for number of coral larvae each point represents

## The dataframe above has one attribute per point: % live coral cover, ranging from 0 - 64.45%.
## Larvae at each point is proportional to coral cover. Since there are no fractional larvae,
## I will first translate coral cover (i.e., larval amount) to the nearest integer.
## Then, for computational efficiency, I will use each larvae to represent ten larvae, so will find the
## nearest integer to 1/10 coral cover.

## Translate coral cover to nearest integer.
sf_40m$larvae <- round(sf_40m$live_coral)
sum(sf_40m$larvae) # There are 174435 larvae in this version

## Same, but if all larvae represent 10 larvae (faster computation)
sf_40m$larv_red <- round(sf_40m$larvae/10)
sum(sf_40m$larv_red) # There are 17051 larvae in this version

## Export to shapefile
st_write(sf_40m, "/Users/rachelcarlson/Documents/Research/Larvae/Data/Primary/gao_pointcloud.shp")

