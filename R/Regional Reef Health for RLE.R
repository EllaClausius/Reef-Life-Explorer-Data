### Calculate Regional Reef Health for RLE 
### EClausius 2020

rm(list = ls())
library(tidyverse)
library(sf)
library(raster)
library(geosphere)


load("data/data-processed/metrics_global.rda")
load("data/data-processed/sdat.rda")


metrics_global <- merge(metrics_global, sdat[, c("SiteCode", "Province")])

IndicatorsSite <- metrics_global %>% 
  group_by(SiteCode, Lat, Lon, Province) %>%
  summarise(B20 = mean(B20), 
            Habitat = mean(habitat), 
            SharksRays = mean(sharks)) %>% 
  ungroup()


#calc quantiles
quantile(IndicatorsSite$SharksRays, probs = c(0.5, 0.75, 0.9, 0.95))
quants <- IndicatorsSite %>%
  summarize_at(vars(B20:SharksRays), quantile, probs = 0.75, na.rm= TRUE)



#new variables that are true/false for if value is greater than quantile 
IndicatorsSite$SharksQ <- IndicatorsSite$SharksRays > quants$SharksRays
IndicatorsSite$HabitatQ <- IndicatorsSite$Habitat > quants$Habitat
IndicatorsSite$B20Q <- IndicatorsSite$B20 > quants$B20



provinces<- IndicatorsSite %>% 
  group_by(Province) %>% 
  summarise(lat = mean(Lat), lon = mean(Lon), nsurveys = n(),
            habitat = mean(Habitat, na.rm=TRUE), B20 = mean(B20), sharks = mean(SharksRays), 
            Q1sharks = 100*sum(SharksQ)/n(),
            Q1B20 = 100*sum(B20Q)/n(),
            Q1Habitat = 100*sum(HabitatQ, na.rm=TRUE)/n())



p1 <- provinces %>% dplyr::select(Province, habitat, B20, sharks) %>% 
  gather(indicator, health_v, -Province)

p2 <- provinces %>% dplyr::select(Province, Q1sharks, Q1B20, Q1Habitat) %>% 
  gather(indicator, health_percentage, -Province) 
p2$indicator[p2$indicator == "Q1sharks"] <- "sharks"
p2$indicator[p2$indicator == "Q1Habitat"] <- "habitat"
p2$indicator[p2$indicator == "Q1B20"] <- "B20"

provinces <- merge(p1, p2, by=c("Province", "indicator"))


load("data/data-processed/indicator_ids.rda")
load("data/data-processed/spatial_levels.rda")

spatial_levels <- spatial_levels %>% subset(type_id == 2)

colnames(provinces)[1] <- "name"

prov_health <- merge(x=provinces, y=spatial_levels[, c("id", "name")], by="name", all.x=TRUE)
prov_health <- merge(prov_health, metric_id, by="indicator", all.x=TRUE)

prov_health[, c(1, 2)] <- NULL
colnames(prov_health)[3] <- "spatial_level_id"

prov_health <- prov_health[, c(3, 4, 1, 2)]
prov_health <- tibble::rowid_to_column(prov_health, "id")
write.csv(prov_health, "outputs/data-outputs/provincial_health.csv", row.names=FALSE)
