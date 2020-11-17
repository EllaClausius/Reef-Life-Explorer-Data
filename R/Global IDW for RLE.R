### Create Global IDW heatmaps for RLE 
## EClausius November 2020 

rm(list = ls())
library(tidyverse)
library(sf)
library(raster)
library(geosphere)


load("data/data-processed/metrics_global.rda")
source("R/Functions for RLE Script.R") ##functions needed for IDW 

select <- dplyr::select

res <- c(1000, 1000) #resolution of map x, y
my_buffer <- 100*1000 # The size of the buffers around the points in metres
p <- 1.6 #power for IDW
      #0.5 was the default in the old code 
      #lower values mean its more influenced by distances further away
      # higher values = less smoothing 

mycolpal <- source("data/data-raw/IDW palette_new.txt")$value  #colour palette for IDW output


#### Set up data input -----------------------------------------------------------------------------------------------------------

global_indic  <- metrics_global %>%
  group_by(SiteCode) %>% 
  summarise(lat = mean(Lat),
            lon = mean(Lon),
            B20 = mean(B20, na.rm = T),
            fish.richness = mean(fish.richness, na.rm=T),
            CF.richness = mean(CF.richness, na.rm=T),
            invert.richness = mean(invert.richness, na.rm=T),
            CTI = mean(CTI, na.rm=T),
            Urchins = mean(Urchins, na.rm=T), 
            sharks = mean(sharks, na.rm=T), 
            habitat = mean(habitat, na.rm=T)) 


#### Set global projection --------------------------------------------------------------------------------------------------------
 
sf_global <- st_as_sf(global_indic, coords = c("lon", "lat"))
st_crs(sf_global) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" 



#### Biodiversity Indicator IDWs  ------------------------------------------------------------------------------------------------


raster.ls <- list()  #set up list of IDW values 

sites <- global_indic %>%  select(SiteCode, lat, lon) %>% na.omit()

sites_xy <- cbind(sites$lon, sites$lat)

## Setup raster and identify cells close to a site (within buffer) ---------------------------------------------------------------
r <- raster(nrow = res[1], ncol = res[2])
icell <- cellFromXY(r, sites_xy)
r[icell] <- 1
rdist <- distance(r)
ifilter <- which(rdist[] < my_buffer)
xyinterp <- xyFromCell(r, ifilter)


xydist <- distm(sites_xy, xyinterp)/1000  #Now get geodesic distance matrix from all included cells to sitesm

#idw_vals <- map(metrics, ~IDW_indicators(.x, xydist, p, global_indic))

metrics <- c("B20", "fish.richness", "CF.richness", "invert.richness", "Urchins", "sharks") 

for(i in 1:length(metrics)){ 
  this_metric <- metrics[i]
  data.sub <- global_indic %>% 
    select(SiteCode, lat, lon, metric = all_of(this_metric)) 
  indicator_values <- pull(data.sub, metric)
  w <- 1/(xydist^p)
  indic <- matrix(indicator_values, nrow = 1)
  isreal <- which(!is.na(indic))
  nvals <- length(isreal)
  
  val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
  val <- as.numeric(val)
  dat_global <- data.frame(xyinterp, val)
  names(dat_global) <- c("lon", "lat", this_metric)
  
  ### assign HEX values -------------------------------------
  probs <-  seq(0, 1, by =0.05)
  
  dathex <- map(this_metric, ~gethex(.x, dat_global, pal = mycolpal, probs = seq(0, 1, by =0.05)))  %>% 
    data.frame()
  names(dathex) <- paste0(this_metric, "_HEX")
  
  dat_global <- bind_cols(dat_global, dathex)
  dat_global$metric <- this_metric
  
  dat_global <- dat_global[, c(2, 1, 5, 3, 4)]
  names(dat_global) <- c("lat", "lon", "indicator", "value", "colour")
  
  raster.ls[[i]] <- dat_global
  
}

global_metrics_raster <- do.call("rbind", raster.ls)



#### Habitat IDW -----------------------------------------------------------------------------------------------------------------------------------------------------------------

raster.ls <- list()
sites <- global_indic %>%  select(SiteCode, lat, lon, habitat) %>% na.omit()
sites_xy <- cbind(sites$lon, sites$lat)
r <- raster(nrow = res[1], ncol = res[2])
icell <- cellFromXY(r, sites_xy)
r[icell] <- 1
rdist <- distance(r)
ifilter <- which(rdist[] < my_buffer)
xyinterp <- xyFromCell(r, ifilter)
xydist <- distm(sites_xy, xyinterp)/1000
metrics <- "habitat"

for(i in 1:length(metrics)){ 
  this_metric <- metrics[i]
  data.sub <- global_indic %>% 
    select(SiteCode, lat, lon, metric = all_of(this_metric)) %>% na.omit()
  indicator_values <- pull(data.sub, metric)
  w <- 1/(xydist^p)
  indic <- matrix(indicator_values, nrow = 1)
  isreal <- which(!is.na(indic))
  nvals <- length(isreal)
  
  val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
  val <- as.numeric(val)
  dat_global <- data.frame(xyinterp, val)
  names(dat_global) <- c("lon", "lat", this_metric)
  
  ### assign HEX values -------------------------------------
  probs <-  seq(0, 1, by =0.05)
  
  dathex <- map(this_metric, ~gethex(.x, dat_global, pal = mycolpal, probs = seq(0, 1, by =0.05)))  %>% 
    data.frame()
  names(dathex) <- paste0(this_metric, "_HEX")
  
  dat_global <- bind_cols(dat_global, dathex)
  dat_global$metric <- this_metric
  
  dat_global <- dat_global[, c(2, 1, 5, 3, 4)]
  names(dat_global) <- c("lat", "lon", "indicator", "value", "colour")
  
  raster.ls[[i]] <- dat_global
  
}

global_habitat_raster <- do.call("rbind", raster.ls)



#### CTI IDW -------------------------------------------------------------------------------------------------------------------------------------------------


raster.ls <- list()
sites <- global_indic %>%  select(SiteCode, lat, lon, CTI) %>% na.omit()
sites_xy <- cbind(sites$lon, sites$lat)
r <- raster(nrow = res[1], ncol = res[2])
icell <- cellFromXY(r, sites_xy)
r[icell] <- 1
rdist <- distance(r)
ifilter <- which(rdist[] < my_buffer)
xyinterp <- xyFromCell(r, ifilter)
xydist <- distm(sites_xy, xyinterp)/1000

metrics <- "CTI"

for(i in 1:length(metrics)){ 
  this_metric <- metrics[i]
  data.sub <- global_indic %>% 
    select(SiteCode, lat, lon, metric = all_of(this_metric)) %>% na.omit()
  indicator_values <- pull(data.sub, metric)
  w <- 1/(xydist^p)
  indic <- matrix(indicator_values, nrow = 1)
  isreal <- which(!is.na(indic))
  nvals <- length(isreal)
  
  val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,])
  val <- as.numeric(val)
  dat_global <- data.frame(xyinterp, val)
  names(dat_global) <- c("lon", "lat", this_metric)
  
  ### assign HEX values -------------------------------------
  probs <-  seq(0, 1, by =0.05)
  
  dathex <- map(this_metric, ~gethex(.x, dat_global, pal = mycolpal, probs = seq(0, 1, by =0.05)))  %>% 
    data.frame()
  names(dathex) <- paste0(this_metric, "_HEX")
  
  dat_global <- bind_cols(dat_global, dathex)
  dat_global$metric <- this_metric
  
  dat_global <- dat_global[, c(2, 1, 5, 3, 4)]
  names(dat_global) <- c("lat", "lon", "indicator", "value", "colour")
  
  raster.ls[[i]] <- dat_global
  
}

global_cti_raster <- do.call("rbind", raster.ls)

#### Combine all ------------------------------------------------------------------------------------------------------------------------


global_metrics_raster_all <- rbind(global_metrics_raster, global_habitat_raster, global_cti_raster)


#### CREATE global_markers.csv table for input to Reef Life Explorer system  -----------------------------------------------------------

load("data/data-processed/indicator_ids.rda")
global_markers <- merge(global_metrics_raster_all, metric_id, "indicator")  ##assign indicator_ids to each row (numeric)
global_markers$indicator <- NULL #then drop indicator column (character)

global_markers <- global_markers[, c(1, 2, 5, 3, 4)]
global_markers <- tibble::rowid_to_column(global_markers, "id")   ##add row Ids (for input into the system)
write.csv(global_markers, "outputs/data-outputs/global_markers.csv", row.names=FALSE)
