library(raster)
library(sf)
library(fasterize)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(here)

loc <- st_read(here("data/GIS/","red_eBird_locations.shp"))

world <- ne_countries(scale = "large", returnclass = "sf")
world <- st_transform(world, crs(loc))

rr <- raster(here("cov_raster/","CR_CL_MN.tif"))
test_spdf <- as(rr, "SpatialPixelsDataFrame")
test_df <- as.data.frame(test_spdf)
colnames(test_df) <- c("value", "x", "y")


gg <- ggplot() +  
  geom_sf(data = world) +
  geom_tile(data=test_df, aes(x=x, y=y, fill=value)) + 
  scale_fill_gradientn(name = "Crown closure [scaled]", colours = c('#d7191c','#fdae61','#ffffbf','#abd9e9','#2c7bb6')) +
  geom_point(data = loc, aes(x = X, y = Y), size = 1, 
             shape = 16, fill = "darkred") +
  xlab("") + ylab("") +
  coord_sf(xlim = c(min(test_df$x), max(test_df$x)), ylim = c(min(test_df$y), max(test_df$y)))
  
ggsave(here("figures/","Point locations plus crown closure.png"), plot = gg, dpi = 300)
