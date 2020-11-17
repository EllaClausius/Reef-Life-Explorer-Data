#### Create Location IDW heatmaps foe RLE 
#### EClausius November 2020 

rm(list = ls())
library(tidyverse)
library(sf)
library(raster)
library(geosphere)
library(rgdal)
library(fasterize)
library(data.table)
select <- dplyr::select

#test
load("data/data-processed/sites-gamms.rda")

size_specs <- dfsites %>% 
  group_by(Location) %>% 
  summarize(max.lat = max(Lat), 
            min.lat = min(Lat), 
            max.lon = max(Lon), 
            min.lon = min(Lon)) %>% 
  group_by(Location) %>% 
  summarize(lat.diff = (max.lat - min.lat), 
            lon.diff = (max.lon - min.lon)) 
size_specs$factor <- ifelse(size_specs$lat.diff > size_specs$lon.diff, "long", "wide")

##source functions we will need
source("R/Functions for RLE Script.R")
load("data/data-processed/loc_sizes.rda")


mycolpal <- source("data/data-raw/IDW palette_new.txt")$value  #colour palette


p <- 1.6 #power for IDW
#0.5 was the default in the old code 
#lower values mean its more influenced by distances further away
# higher values = less smoothing 

probs <- seq(0, 1, by =0.05) #quantile intervals for hex colours


aus <- readOGR(dsn = "data/data-raw/shapefiles/aus-shapefile", "AUS_adm1") #%>% spTransform(raster::crs(local_crs))#bring in Aus Shapefile


metrics <- c("Urchins", "fish.richness", "B20", "CTI", "invert.richness", "CF.richness") 

dfout <- NULL
rout <- NULL
routnames <- NULL

res <- c(100, 150, 200, 400, 800)
size <- c("X-small", "Small", "Medium", "Large", "X-large")
buffers <- c(500, 1000, 2000, 4000, 8000)
size.buffers <- data.frame(size, buffers,res)

## RUN FOR TESTING 

#this_size <- "Small"
#locations <- "Lord Howe Island"
#my_buffer <- 1000
#metrics <- "CTI"

##
mapbox_testing <- FALSE    ##spits out csvs of each local IDW for testing and zoom/size determination 



####  Run IDW: loop over size classes, then locations, then metrics then years. -----------------------------------------------------


##Loop over size classes 
for(i_size in 1:length(size)){
  this_size <- size[i_size]
  #this_size <- "Medium"  ##for testing
  size.locs <- loc_sizes$Location[loc_sizes$size == this_size]
  locations <- unique(dfsites$Location)
  locations <- locations[locations %in% size.locs]
  my_buffer <- size.buffers$buffers[size == this_size]
  #locations <- c("Port Phillip Bay", "Port Phillip Heads", "Rottnest Island", "Sydney", "Bicheno")
  #print(locations)
  
  #Loop over locations 
  for (i_location in 1:length(locations)){
    this_location <- locations[i_location]
    #this_location <- locations[1]
    print(this_location)
    
    dftemp <- filter(dfsites, (Location == this_location))
    years <- unique(dftemp$year_month)
    
    #
    #Location definition and distance matrix 
    #
    dfloc <- dftemp %>% dplyr::select(SiteCode, Lon, Lat) %>% distinct()
    sites_xy <- cbind(dfloc$Lon, dfloc$Lat)
    
    
    lon.min <- min(dftemp$Lon)
    lon.max <- max(dftemp$Lon)
    lat.min <- min(dftemp$Lat)
    lat.max <- max(dftemp$Lat)
    lon.diff <- lon.max - lon.min
    lat.diff <- lat.max - lat.min
    diff.wide <- lat.diff/lon.diff
    diff.long <- lon.diff/lat.diff
    cols <- ifelse(lon.diff > lat.diff, size.buffers$res[size.buffers$size == this_size], 
                   (size.buffers$res[size.buffers$size == this_size])*diff.long) 
    
    rows <- ifelse(lon.diff > lat.diff, (size.buffers$res[size.buffers$size == this_size])*diff.wide, 
                   size.buffers$res[size.buffers$size == this_size]) 
    
    
    
    r <- raster(extent(c(lon.min, lon.max,
                         lat.min, lat.max)),
                nrow = rows, ncol = cols)
    
    proj4string(r) <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0" 
    r <- extend(r, ceiling((res(r)*110*1000)/my_buffer)*150) ##extends buffer around local extent. Increase last value if IDW cutting off edges. 
    
    
    icell <- cellFromXY(r, sites_xy) ### gets r cells that match coordinates from sites_xy
    r[icell] <- 1 ###assigns grid cells from icell a value of 1 
    rdist <- distance(r) ###computes the distance, for all cells that are NA, to the nearest cell that is not NA (i.e., to the sites in sites_xy)
    #save(rdist, file = 'data/rdist.rda', overwrite = T)
    ifilter <- which(rdist[] < my_buffer)  # a list of the grid cells in rdist that fall within the defined buffer 
    xyinterp <- xyFromCell(r, ifilter) ##generate a list of coordinates for each of the grid cells that fall within the buffer of each site
    xydist <- distm(sites_xy, xyinterp)/1000 #calculate the distance of each grid cell inside the buffer from each site
    w <- 1/(xydist^p) #weights each grid cell value by it's distance from the site according to the power you assign. 
    
    #
    # IDW calculation by each metric and each year 
    #
    #Loop over metrics
    for (i_metrics in 1:length(metrics)){
      this_metric <- metrics[i_metrics]
      #this_metric <- "CTI"
      print(this_metric)
      
      #Loop over years
      for (i_years in 1:length(years)){
        this_year <- years[i_years]
        
        dfyear <- filter(dftemp, (year_month == this_year) & (metric == this_metric))
        dfyear <- na.omit(dfyear)
        
        indicator_values <- pull(dfyear, "fit")  ##pulls all indicator values from dftemp
        indic <- matrix(indicator_values, nrow = 1) #creates a matrix with 1 row containing indicator values
        isreal <- which(!is.na(indic)) #determines which matric columns are not NA 
        #nvals <- length(isreal) #number of real values 
        val <- (indic[isreal] %*% w[isreal,])/colSums(w[isreal,]) ##multiples the distance matrix cells by the real indicator values to get the interpolated values. 
        
        rval <- r
        rval[ifilter] <- as.numeric(val)
        
        this_raster <- rval 
        aus1 <- crop(aus, extent(this_raster)) #crops the raster to the extent of the polygon to improve the speed of the mask 
        
        this_raster2 <- raster::mask(this_raster, aus1, inverse=TRUE) #Removes the part of the this_raster that overlaps with the cropped aus1 rasterr)
        raster.df <- rasterToPoints(this_raster2)  #convert the raster to a dataframe 
        dfvalues <- data.frame(Location = this_location, 
                               Metric = this_metric, 
                               Year = this_year, 
                               Lon = raster.df[,1],
                               Lat = raster.df[, 2], 
                               prediction = as.numeric(raster.df[, 3]))
        
        dfout <- c(dfout, list(dfvalues))
      }
    }
  }
}

###if this spits out an "Error in validityMethod(object): invalid extent: xmin >= xmax" error, it's because there are no 
###Large locations in the dfsites dataframe. 


dfx <- do.call("rbind", dfout)  




#### Assign hex values to each point ------------------------------------------------------------------------------------------------


gethexlocal <- function(x, probs, pal = mycolpal){
  xcols <- leaflet::colorQuantile(pal= mycolpal, 
                                  x, reverse = TRUE,
                                  probs = probs)
  y <- try(xcols(x), 
           TRUE)
  if(class(y) != "try-error"){
    y
  } else {
    mycolpal[21]
  }
}

locations2 <- unique(dfx$Location)
dfx$hex <- NA
for (i_location in 1:length(locations2)){
  this_location <- locations2[i_location]
  for (i_metrics in 1:length(metrics)){
    this_metric <- metrics[i_metrics]
    irow <- which((dfx$Location == this_location) & 
                    (dfx$Metric == this_metric))
    dfx$hex[irow] <- gethexlocal(dfx[irow,]$prediction, probs, mycolpal)
  }
}



#### Check by plotting  -----------------------------------------------------------------------------------------------------------


idw_plots <- FALSE

locations <- unique(dfx$Location)
if(idw_plots){
  for(i_location in 1:length(locations)){
    this_location <- locations[i_location]
    data <- dfx[dfx$Location == this_location,]
    max_year <- max(data$Year)
    data <- data[data$Year == max_year,]
    p <- ggplot(data) + 
      geom_point(aes(Lon, Lat), colour=data$hex, size=0.5)
    ggsave(p, width = 10, height = 10,
           filename = paste0("idw-plots/", this_location, ".png"))
  }
} else{
  print("No idw-plot outputs")
}

#### Year testing 
#ppb <- filter(dfx, Location == "" & Metric == "CTI")
#ggplot() + geom_point(data=ppb, aes(x=Lon, y=Lat, color=prediction)) + facet_wrap(~Year)


#### Create datatables for mapbox testing -----------------------------------------------------------------------------------------


locations <- unique(dfsites$Location)

if(mapbox_testing){ 
  
  for (i_location in 1:length(locations)){
    this_location <- locations[i_location]
    
    sub.dfx <- dfx[dfx$Location == this_location,]
    sub.dfx <- sub.dfx[sub.dfx$Metric == "B20",]
    max_year <- max(sub.dfx$Year)
    sub.dfx <- sub.dfx[sub.dfx$Year == max_year, ]
    sub.dfx[, c(1:3)] <- NULL 
    write.csv(sub.dfx, file=paste("local-IDW-testing/",this_location, ".csv", sep=""), row.names = FALSE) 
  }
} else{
  print("No mapbox testing")
}


if(mapbox_testing){ 
  
  for (i_size in 1:length(size)){
    this_size <- size[i_size]
    loc.this_size <- loc.sizes[loc.sizes$size == this_size, ]
    
    sub.dfx <- dfx[dfx$Location %in% loc.this_size$Location,]
    sub.dfx <- sub.dfx[sub.dfx$Metric == "B20",]
    
    #max_year <- max(sub.dfx$Year)
    sub.dfx <- sub.dfx[sub.dfx$Year == 2018, ]
    sub.dfx[, c(2:3)] <- NULL 
    write.csv(sub.dfx, file=paste("local-IDW-testing/",this_size, ".csv", sep=""), row.names = FALSE) 
  }
} else{
  print("No mapbox testing")
}







#### CREATE location_markers.csv table for input into Reef Health Explorer system ---------------------------------------------------

l_markers <- dfx

load("data/data-processed/indicator_ids.rda")
load("data/data-processed/spatial_levels.rda")

colnames(l_markers) <- c("name", "indicator", "year", "lon", "lat", "value", "color")


l_markers2 <- merge(l_markers, metric_id, "indicator")
spatial_levels <- spatial_levels %>%  subset(type_id == 4)
l_markers2 <- merge(x=l_markers2, y=spatial_levels[, c("id", "name")], by="name", all.x=TRUE)
l_markers2[, c(1, 2)] <- NULL
colnames(l_markers2)[7] <- "spatial_level_id"

location_markers <- l_markers2[, c(3, 2, 1, 7, 6, 4, 5)]
location_markers <- tibble::rowid_to_column(location_markers, "id")
write.csv(location_markers, "outputs/data-outputs/location_markers.csv", row.names=FALSE, quote=F)
# write.table(location_markers, location_markers, sep=",", col.names=FALSE)


#### Subset location_markers by indicator for easier data import --------------------------------------------------------------------

location_markers_B20 <- subset(location_markers, indicator_id == 4)
write.csv(location_markers_B20, "outputs/data-outputs/location-markers/location_markers_B20.csv", row.names=FALSE, quote=F)

location_markers_fish <- subset(location_markers, indicator_id == 1)
write.csv(location_markers_fish, "outputs/data-outputs/location-markers/location_markers_fish.csv", row.names=FALSE, quote=F)

location_markers_CF <- subset(location_markers, indicator_id == 6)
write.csv(location_markers_CF, "outputs/data-outputs/location-markers/location_markers_CF.csv", row.names=FALSE, quote=F)

location_markers_invert <- subset(location_markers, indicator_id == 2)
write.csv(location_markers_invert, "outputs/data-outputs/location-markers/location_markers_invert.csv", row.names=FALSE, quote=F)

location_markers_urchins <- subset(location_markers, indicator_id == 3)
write.csv(location_markers_urchins, "outputs/data-outputs/location-markers/location_markers_urchins.csv", row.names=FALSE, quote=F)

location_markers_CTI <- subset(location_markers, indicator_id == 5)
write.csv(location_markers_CTI, "outputs/data-outputs/location-markers/location_markers_CTI.csv", row.names=FALSE, quote=F)



