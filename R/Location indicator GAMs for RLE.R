### Run location-level GAMs for RLE 
## EClausius November 2020 

rm(list = ls())
library(tidyverse)
library(mgcv)
library(dplyr)
library(cowplot)

select <- dplyr::select
load("data/data-processed/metrics_global.rda")

#Min years and sites for a location to be included
minyears <- 3
minsites <- 5
maxyear <- 2015 #only include locations that have been surveyed since at least 2015
resurvey_rule <- TRUE #do you want to require the same minsites sites having to be resurveyed minyears times?

save_plots <- FALSE #Do you want to save plots (will do so in subfolders of the current wd)


metrics <- c("B20", "fish.richness", "CF.richness", "invert.richness", "CTI", "Urchins")

source("R/Functions for RLE Script.R")
metrics_global <- metrics_global %>% 
  mutate(Year = lubridate::year(survey_date))


#### Identify locations that meet data requirements ------------------------------------------------------------------------------

if(resurvey_rule){
  #if sites do need to be resurveyed
  location_include <- metrics_global %>% 
    select(SiteCode, `Country`,Location, Year) %>%
    distinct() %>%
    group_by(Location, `Country`, SiteCode) %>%
    summarize(nyears = n_distinct(Year), 
              lastyear = max(Year)) %>%
    group_by(Location) %>% 
    summarize(max_resurvey = max(nyears),
              nsites = n_distinct(SiteCode), 
              most_recent_year = max(lastyear)) %>%
    filter((max_resurvey > minyears) & (nsites > minsites)) %>% 
    filter(most_recent_year >= maxyear)
  
} else{
  #If sites don't have to be resureyed
  met_sum <- metrics_global %>% 
    select(SiteCode, `Country`,Location, Year) %>%
    distinct() %>%
    group_by(Location, `Country`) %>%
    summarize(nyears = n_distinct(Year), 
              nsites = n_distinct(SiteCode))
  location_include <- filter(met_sum, (nyears >= minyears) & (nsites >= minsites))
}


#
# Filter locations in metrics
#

#Make sure to do this by country AND location, because some location names are duplicated across countries 

logplus <- min(metrics_global$B20[metrics_global$B20>0], na.rm = TRUE)/10

metrics_location <- metrics_global %>%
  inner_join(location_include) %>%
  mutate(month = lubridate::month(survey_date),
         year_month = Year + month/12) %>%
  data.frame()

# table(metrics_location$`Country/Region`)
# table(metrics_location$Location)


# ----------------- #
# GAMs
# ----------------- #

#
#Fit Location by Location as initial test
#

all_locations <- metrics_location %>%  pull(Location) %>% unique()  #RUN ALL LOCATIONS


for (i in 1:length(metrics)) dir.create(paste0("outputs/plots-gams/",metrics[i]))#create directories for storing graphs in


# Loop over all metrics (map2 was causing me issues)
xout <- NULL
for (i in 1:length(metrics)){
  print(metrics[i])
  xtemp <- map(all_locations, 
               ~fitgam(.x, metrics_location, metrics[i], 
                       save_plots = save_plots))
  xout <- c(xout, list(xtemp))
}

#
#Export mean GAM fits 
#
x <- flatten(xout)
errorlist <- map(x, class) %>% unlist()
z <- NULL
for (i in which(errorlist != "try-error")) z <- c(z, list(x[[i]]$predictions))

df <- do.call("rbind", z) %>% 
  dplyr::select(-SiteCode) %>%
  dplyr::select(Location, Metric = metric, Year = year_month, 
                Mean = fit, StandardError = se)


#### Export site by site fits  ---------------------------------------------------------------------------------------------------

zsites <- NULL
for (i in which(errorlist != "try-error")) zsites <- c(zsites, list(x[[i]]$site_predictions))

dfsites <- do.call("rbind", zsites) 


#### Location data - site means --------------------------------------------------------------------------------------------------

site_year_means <- metrics_location %>% 
  group_by(SiteCode, Year, Country) %>%
  summarize_at(metrics, mean, na.rm = TRUE)


#### Save outputs  -----------------------------------------------------------------------------------------------------------------


#
# #output for web interface
# write_csv(df, "data/Table_4-GAMM-output-locations.csv")
# write_csv(site_year_means, "data/Table_5-Site-by-year-means.csv")

##output for IDW calculations 
save(dfsites, file = "data/data-processed/sites-gamms.rda")


#### CREATE location_trend.csv table for input into Reef Life Explorer system  ------------------------------------------------------


load("data/data-processed/indicator_ids.rda")
load("data/data-processed/spatial_levels.rda")

l_trend <- df

#\ l_trend$Mean <- ifelse(l_trend$Metric == "CTI", as.numeric(l_trend$Metric, ))
# l_trend$Mean <- ifelse(l_trend$Mean[l_trend$Metric %in% c("fish.richness", "CF.richness", "invert.richness", "Urchins", "B20")] <= 0, 0, as.numeric(l_trend$Metric))
# l_trend$Mean <- case_when(l_trend$Metric == "CTI" ~ as.numeric(Mean), 
#                           l_trend$Mean <= 0 ~ 0, 
#                           TRUE ~ as.numeric(l_trend$Mean))

l_trend$min_v <- l_trend$Mean - l_trend$StandardError
l_trend$max_v <- l_trend$Mean + l_trend$StandardError
l_trend$StandardError <- NULL 

colnames(l_trend)[c(2, 4, 1)] <- c("indicator", "mean_v", "name")

l_trend2 <- merge(l_trend, metric_id, "indicator")
spatial_levels <- spatial_levels %>%  subset(type_id == 4)
l_trend2 <- merge(x=l_trend2, y=spatial_levels[, c("id", "name")], by="name", all.x=TRUE)
l_trend2$name <- NULL
l_trend2$indicator <- NULL
colnames(l_trend2)[6] <- "spatial_level_id"
colnames(l_trend2)[1] <- "year"


location_trend <- l_trend2[, c(3, 2, 4, 1, 6, 5)]
location_trend <- na.omit(location_trend)
location_trend <- location_trend[!(location_trend$indicator_id == "9"), ]
location_trend <- tibble::rowid_to_column(location_trend, "id")
write.csv(location_trend, "outputs/data-outputs/location_trend.csv", row.names=FALSE)

