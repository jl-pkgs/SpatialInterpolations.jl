pacman::p_load(
  Ipaper, data.table, dplyr, lubridate, 
  terra, sf2
)

ra = rast("./Project_十堰/data/dem_ShiYan_1km.tif")
ra2 = aggregate(ra, fact=2)

writeRaster(ra2, "./Project_十堰/data/dem_ShiYan_2km.tif")
