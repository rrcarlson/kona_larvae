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
points_2m <- rasterToPoints(gao)



