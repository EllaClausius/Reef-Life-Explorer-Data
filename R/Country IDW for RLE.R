### Create Country IDW for RLE 
## EClausius November 2020

rm(list = ls())
library(tidyverse)
library(sf)
library(raster)
library(geosphere)


load("data/data-processed/metrics_global.rda")

#source functions we will need
source("R/Functions for RLE Script.R")

#res <- c(800, 800) #resolution of map x, y
my_buffer <- 100000 # The size of the buffers around the points in metres
p <- 2 #power for IDW
#0.5 was the default in the old code 
#lower values mean its more influenced by distances further away
# higher values = less smoothing 

#Colour palette to use for outputing hex codes. 
mycolpal <- source("data/data-raw/IDW palette.txt")$value

#
###Metric variable names to make summaries for
#
metrics <- c("B20", "fish.richness", "CF.richness", "invert.richness", "CTI", "Urchins", "sharks")



### Average across time 

colnames(metrics_global)[6] <-  "Country"

aus_indic  <- metrics_global %>% 
  filter(Country == "Australia") %>%  ##change this when more countries become available 
  group_by(SiteCode) %>% 
  # mutate(B20_log = log10(B20+0.001)) %>% 
  summarise(lat = mean(Lat),
            lon = mean(Lon),
            B20 = mean(B20, na.rm = T),
            #log10B20 = mean(log10B20, na.rm = T),
            fish.richness = mean(fish.richness, na.rm=T),
            CF.richness = mean(CF.richness, na.rm=T),
            invert.richness = mean(invert.richness, na.rm=T),
            CTI = mean(CTI, na.rm=T),
            Urchins = mean(Urchins, na.rm=T), 
            sharks = mean(sharks, na.rm=T))

sf_aus <- st_as_sf(aus_indic, coords = c("lon", "lat"))
st_crs(sf_aus) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" 


#
# IDW raster grid and distances
#
sites_xy <- cbind(aus_indic$lon, aus_indic$lat)


lon.min <- min(aus_indic$lon)
lon.max <- max(aus_indic$lon)
lat.min <- min(aus_indic$lat)
lat.max <- max(aus_indic$lat)
lon.diff <- lon.max - lon.min
lat.diff <- lat.max - lat.min
diff <- lat.diff/lon.diff


#setup raster and identify cells close to 
# a site (within buffer)
r <- raster(nrow = 1500*diff, ncol = 1500)
icell <- cellFromXY(r, sites_xy)
r[icell] <- 1
rdist <- distance(r)
#
ifilter <- which(rdist[] < my_buffer)
xyinterp <- xyFromCell(r, ifilter)

#Now get geodesic distance 
# matrix from all included cells to sites
xydist <- distm(sites_xy,
                xyinterp)/1000


idw_vals <- map(metrics, ~IDW_indicators(.x, xydist, p, aus_indic))

dat_aus <- data.frame(xyinterp, do.call("cbind", idw_vals))
names(dat_aus) <- c("lon", "lat", metrics)
dat_aus2 <- dat_aus

# rval <- r
# rval[ifilter] <- dat_global$B20
# plot(rval)

#
### Colour scale 
#

dathex <- map(metrics, ~gethex(.x, dat_aus, pal = mycolpal, probs = seq(0, 1, by =0.05))) %>%
  do.call("cbind", .) %>% data.frame()


names(dathex) <- paste0(metrics, "_HEX")
dathex <- dathex %>%  gather(indicator, colour, 1:7)
dat_aus <- dat_aus %>%  gather(indicator, value, -lat, -lon)
dat_aus <- bind_cols(dat_aus, dathex)
dat_aus2 <- dat_aus[,c(1, 2, 3, 4, 6) ]

##################################################


#
### Check by plotting 
#


# par(bg = "#192638")
# plot(dat_aus$lon, dat_aus$lat, pch = 15,  cex=0.1, col = as.character(dat_aus$B20_HEX))
# plot(dat_aus$lon, dat_aus$lat, pch = 15, cex = 0.4, col = as.character(dat_aus$fish.richness_HEX))
# plot(dat_aus$lon, dat_aus$lat, pch = 15, col = as.character(dat_aus$CF.richness_HEX))
# plot(dat_aus$lon, dat_aus$lat, pch = 15, cex = 1, col = as.character(dat_aus$CTI_HEX))
# plot(dat_aus$lon, dat_aus$lat, pch = 15, col = as.character(dat_aus$Urchins_HEX))
# plot(dat_aus$lon, dat_aus$lat, pch = 15, col = as.character(dat_aus$habitat_HEX))


#
### CREATE global_markers.csv table for input into the Reef Life Explorer system 
#

aus_markers <- dat_aus2
load("data/data-processed/indicator_ids.rda")
aus_markers <- merge(aus_markers, metric_id, "indicator")

aus_markers$indicator <- NULL

aus_markers['spatial_level_id']='50'
aus_markers["year"]="2020"
aus_markers <- aus_markers[, c(2, 1, 7, 6, 5, 4)]
aus_markers <- tibble::rowid_to_column(aus_markers, "id")


write.csv(aus_markers, "outputs/data-outputs/country_markers.csv", row.names=FALSE)
