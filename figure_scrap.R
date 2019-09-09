library(raster)
library(sf)
library(fasterize)
library(ggplot2)

setwd("D:/Work_unchanged/NPLCC/eBird/Peter_Arcese_bundle/GIS")

loc <- st_read("red_eBird_locations.shp")

setwd("D:/Work_unchanged/NPLCC/Covariates/NPLCC/GIS/")
covs <- st_read("NPLCC_region_final_100m_grid_covars_scaled.shp")


setwd("D:/Work/Cornell/STEM_downscale/Tradeoffs-biodiversity-cost/figures/")


rr <- raster(covs, res=100)

rr <- fasterize(covs, rr, field = "CR_CL_MN", fun = "max")

test_spdf <- as(rr, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")


gg <- ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value)) + 
  scale_fill_gradientn(name = "Crown closure [scaled]", colours = c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')) +
  geom_sf() +
  geom_point(data = loc, aes(x = X, y = Y), size = 1, 
             shape = 16, fill = "darkred") +
  xlab("") + ylab("")


ggsave("Point locations plus crown closure.png", plot = gg, dpi = 300)






ggplot() +
  geom_sf() +
  geom_sf(data = covs, aes(fill = CR_CL_MN)) +
  geom_sf(data = loc, size = 1, shape = 16, fill = "darkred")   


plot(covs[["CR_CL_MN"]])
plot(loc$geometry, add=TRUE)


test_spdf <- as(rr, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")


ggplot() +  
  geom_tile(data=test_df, aes(x=x, y=y, fill=value), alpha=0.8) + 
  
ggplot()  +
  geom_sf() +
  geom_point(data = loc, aes(x = X, y = Y), size = 4, 
             shape = 23, fill = "darkred")

  geom_polygon(data=OR, aes(x=long, y=lat, group=group), 
               fill=NA, color="grey50", size=0.25) +
  scale_fill_viridis() +
  coord_equal() +
  theme_map() +
  theme(legend.position="bottom") +
  theme(legend.key.width=unit(2, "cm"))
  
  
  