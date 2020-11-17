##Creation of Dataframes for the Reef Life Explorer
###EClausius November 2020 

rm(list = ls())
library(readr)
library(tidyverse)
library(data.table)
library(dplyr)

select <- dplyr::select

load("data/data-processed/fdat.rda")

traits <- fread("data/data-raw/TRAITS MASTER (last update 070720).csv") %>%  as_tibble()
names(traits)[2] <- "taxon"
save(traits, file="data/data-processed/traits.rda")

#### Create Site-level data -----------------------------------------------------------------------------------------------

locdata <- fread("data/data-raw/locations_by country_by province.csv") %>%  as_tibble()
sdat <- fread("data/data-raw/site-masters/Sites-master_consolidated 05112020.csv") %>% as_tibble() %>% 
  select(SiteCode = SiteCode_rev, Latitude = SiteLatitude_rev, Longitude = SiteLongitude_rev, Location = Loc_new) %>% 
  distinct() %>% 
  left_join(locdata, by="Location", all=TRUE) %>%    
  select(SiteCode, Latitude, Longitude, Location = Location_RLE, Country, Province)

sdat <- sdat[!duplicated(sdat$SiteCode),] ##drop SiteCode duplications 
nrow(sdat) == length(unique(sdat$SiteCode)) ##check there are no SiteCode duplications 

save(sdat, file = 'data/data-processed/sdat.rda', overwrite = T)



#### Read in raw RLS M1, M2(fish) and M2(invert) extracts ------------------------------------------------------------=====

m1_in <- fread("data/data-raw/rls-extracts/july2020/ep_m1.csv", na = c("", "NULL")) %>%  as_tibble() ##add updated M1 extracts here 
m2_fish_in <- fread("data/data-raw/rls-extracts/july2020/ep_m1.csv", na = c("", "NULL")) %>%  as_tibble()##add updated M2 Cryptic fish extracts here 
m2_invert_in <- fread("data/data-raw/rls-extracts/july2020/ep_m2_inverts.csv", na = c("", "NULL")) %>%  as_tibble()##add updated M2 Invert extracts here 


m1_in <- subset(m1_in, method != "10") ##get rid of method 10 ATRC Seagrass surveys  

m2_invert_in$biomass <- "NA"  #add biomass column with NAs 
m2_invert_in$method <- 2  #add method column 
m2_fish_in$method <- 2  #add method column

rls_in <- do.call("rbind", list(m1_in, m2_fish_in, m2_invert_in))


    ### Remove duplicated data -------------------------------------------------------------------------------------------------

vic_dups <- fread("data/data-raw/vic_duplications.csv") %>%  as_tibble() 

wa_atrc <- rls_in %>% 
  filter(rls_in$program == "ATRC") %>% 
  filter(area == "Western Australia") %>% 
  dplyr::select(survey_id)###list of all western australian ATRC surveys 

atrc_dups <- rls_in %>% 
  mutate(Year = lubridate::year(survey_date)) %>% 
  filter(rls_in$program == "ATRC") %>% 
  filter(location == "Wilsons Promontory") %>% 
  filter(Year %in% c(1999, 2000)) %>% 
  dplyr::select(!Year) %>% 
  dplyr::select(survey_id)##subset ATRC Wilsons Prom data from 1999 and 2000,then drop it from rls_in due to duplications with Parks Vic data. 

drop.surveys <- rbind(wa_atrc, atrc_dups, vic_dups)
drop.surveys <- distinct(drop.surveys)

rls_in <- rls_in[!(rls_in$survey_id %in% drop.surveys$survey_id), ]  ##remove all western Australian ATRC surveys as requested and duplicated Wilsons Prom data from 1999/2000.

rls_in$biomass <- as.numeric(rls_in$biomass)
rls_in$biomass <- rls_in$biomass/1000 #convert biomass to kg

rls_in <- rls_in %>% 
  select(survey_id, site_code, latitude, longitude, survey_date, program, Location = location, 
         method, block, phylum, class, family, recorded_species_name, species_name, taxon, reporting_name, 
         size_class, total, biomass) %>% 
  distinct() %>% 
  left_join(locdata, by="Location") %>% 
  select(!Location)

names(rls_in)[19] <- "Location"


save(rls_in, file="data/data-processed/rls_in.rda", overwrite=T)


#### Calculate B20 Manually  ------------------------------------------------------------------------------------------------------

B20 <- rls_in %>%   
  filter(method == 1 & class %in% c("Actinopterygii", "Elasmobranchii")) %>% 
  filter(size_class >= 20) %>%  
  group_by(survey_id) %>% 
  summarize(B20 = sum(biomass, na.rm=TRUE)) %>% 
  select(SurveyID = survey_id, B20) 


#### Calculate fish richness manually  --------------------------------------------------------------------------------------------


fish.richness <- rls_in %>% 
  filter(method == 1 & class %in% c("Actinopterygii", "Elasmobranchii", "Chondrichthyes")) %>% 
  group_by(survey_id) %>% 
  summarise(fish.richness = uniqueN(species_name)) %>% 
  group_by(survey_id) %>% 
  summarise(fish.richness = mean(fish.richness)) %>% 
  select(SurveyID = survey_id, fish.richness)


#### Calculate CF.richness manually -----------------------------------------------------------------------------------------------


CF.richness <- rls_in %>% 
  filter(method == 2) %>% 
  filter(class %in% c("Actinopterygii", "Elasmobranchii", "Chondrichthyes")) %>% 
  group_by(survey_id, block) %>% 
  summarise(CF.richness = uniqueN(species_name)) %>% 
  group_by(survey_id) %>% 
  summarise(CF.richness = mean(CF.richness)) %>% 
  select(SurveyID = survey_id, CF.richness)


#### Calculate invert.richness manually ---------------------------------------------------------------------------------------------


invert.richness <- rls_in %>% 
  filter(method == 2 & class %in% c("Asteroidea", "Bivalvia","Cephalopoda",
                                    "Crinoidea", "Echinoidea","Gastropoda", "Holothuroidea", "Malacostraca")) %>% 
  group_by(survey_id, block) %>% 
  summarise(invert.richness = uniqueN(species_name)) %>% 
  group_by(survey_id) %>% 
  summarise(invert.richness = mean(invert.richness)) %>% 
  select(SurveyID = survey_id, invert.richness)


#### Calculate Urchin density manually ----------------------------------------------------------------------------------------------


.getgenus <- function(x){
  unlist(lapply(strsplit(x, split = "[ ]"), function(y) y[1]))
}

urchins <- rls_in %>% 
  filter(method == 2 & class == "Echinoidea") %>%
  mutate(genus = .getgenus(species_name)) %>% 
  filter(genus != "Echinostrephus") %>% 
  group_by(survey_id, block) %>% 
  summarise(urchins = sum(total))  %>% 
  group_by(survey_id) %>% 
  summarise(Urchins = mean(urchins)) %>% 
  select(SurveyID = survey_id, Urchins)



#### Calculate CTI manually ---------------------------------------------------------------------------------------------------------

cti <- rls_in %>% 
  inner_join(traits[, c(2, 18, 27)], by = "taxon") %>% 
  filter(method == 1 & class %in% c("Actinopterygii", "Elasmobranchii", "Chondrichthyes")) %>% 
  # filter(`Water column` != "pelagic non-site attached") %>%  
  mutate(log10Abundance = log10(total + 1), 
         mid.pointXn = log10Abundance*ThermalMP_5min_95max) %>% 
  na.omit() %>% 
  group_by(survey_id) %>% 
  summarise(x = sum(log10Abundance), y = sum(mid.pointXn))  %>% 
  filter(x != 0) %>% 
  na.omit() %>% 
  group_by(survey_id) %>% 
  summarise(CTI = y/x) %>% distinct()  %>% 
  select(SurveyID = survey_id, CTI)



#### Calculate Shark/Ray density  --------------------------------------------------------------------------------------------------

sharks <- rls_in %>% 
  filter(method == 1 & class %in% c("Elasmobranchii", "Chondrichthyes")) %>% 
  group_by(survey_id) %>% 
  summarise(sharks = sum(total, na.rm=FALSE)) %>% 
  select(SurveyID = survey_id, sharks)



#### Calculate Habitat % ---------------------------------------------------------------------------------------------------------

hab_in <- fread("data/data-raw/habitat-data/2020-09-09_Combined extract.csv") %>%  as_tibble() %>% 
          filter(SurveyID != 0) 
names(hab_in)[19] <- "GE_Major_Category"

save(hab_in, file="data/data-processed/hab_in.rda")

hab_totals <- hab_in %>%  
  select(SurveyID, GE_Major_Category, Count) %>% 
  group_by(SurveyID) %>% 
  summarise(total.count = sum(Count))

hab <- hab_in %>%  
  select(SurveyID, GE_Major_Category, Count) %>% 
  filter(GE_Major_Category %in% c("Coral", "Kelp")) %>% 
  group_by(SurveyID, GE_Major_Category) %>% 
  summarise(combined_cover = sum(Count, na.rm=TRUE)) %>% 
  left_join(hab_totals, by="SurveyID", all=T) %>% 
  group_by(SurveyID, GE_Major_Category) %>% 
  summarize(total.cover = (combined_cover/total.count)*100) %>% 
  pivot_wider(names_from = GE_Major_Category, values_from = total.cover)

hab$Kelp <- ifelse(is.na(hab$Kelp), 0, as.numeric(hab$Kelp))
hab$Coral <- ifelse(is.na(hab$Coral), 0, as.numeric(hab$Coral))

hab_in <- hab %>% 
  mutate(habitat = sum(Kelp + Coral))



#### Create Survey-level data (survdat)  -----------------------------------------------------------------------------------------

survdat <- rls_in %>% 
  select(SiteCode = site_code, SurveyID = survey_id, SurveyDate = survey_date) %>% 
  distinct()

survdat <- survdat[!duplicated(survdat$SurveyID),] ##drop SiteCode duplications 
nrow(survdat) == length(unique(survdat$SurveyID)) ##check there are no SiteCode duplications 

save(survdat, file = 'data/data-processed/survdat.rda', overwrite = T) #write survdat rdata file


#### Organise additional bits of data ----------------------------------------------------------------------------------------------


metric_id <- fread("data/data-raw/indicator_ids.csv") %>%  as_tibble()#list of indicator IDs
save(metric_id, file="data/data-processed/indicator_ids.rda", overwrite=T)


spatial_levels <- fread("data/data-raw/spatial-levels/spatial_levels_MASTER 15-09-2020.csv") %>%  as_tibble() #spatial levels for assigning spatial_level_IDs later on
save(spatial_levels, file="data/data-processed/spatial_levels.rda", overwrite=T)


loc_type <- fread("data/data-raw/location_type.csv") %>%  as_tibble() #Location types based on spatial_levels.csv 
save(loc_type, file="data/data-processed/location_type.rda", overwrite=T)

loc_sizes <- fread("data/data-raw/location_sizes.csv") %>%  as_tibble() 
save(loc_sizes, file="data/data-processed/loc_sizes.rda")


#### Creates final surveys.csv table for input into Reef Life Explorer system --------------------------------------------------------

s1 <- survdat ##creates new object with SiteCode, SurveyID & SurveyDate 
colnames(s1) <- c("site_code", "survey_id", "surveydate") 

s2 <- s1 %>%  mutate(Year = lubridate::year(surveydate)) 
s3 <- merge(x=s2, y=spatial_levels[, c("id", "site_code")], by="site_code", all.x=TRUE) ##add spatial_levels id for each SiteCode

surveys <- s3[, c(2, 4, 5)] #reorder columns to fit server requirements 
colnames(surveys)[3] <- "spatial_level_id"

surveys <- surveys[!duplicated(surveys$survey_id),] ##drop survey_id duplications 
nrow(surveys) == length(unique(surveys$survey_id)) ##check there are no SiteCode duplications 
surveys <- na.omit(surveys)  

write.csv(surveys, "outputs/data-outputs/surveys.csv", row.names=F)



#### Create metrics_global table ---------------------------------------------------------------------------------------------------

names(surveys)[1] <- "SurveyID"
survdat2 <- rls_in %>%  
  select(SurveyID = survey_id, Country, SiteCode = site_code, Lat = latitude, Lon = longitude, survey_date, Location) %>% 
  distinct()

metrics_global <- surveys %>%  
  left_join(survdat2,by="SurveyID") %>% 
  full_join(B20, by="SurveyID") %>% 
  full_join(fish.richness, by="SurveyID") %>% 
  full_join(CF.richness, by="SurveyID") %>% 
  full_join(invert.richness, by="SurveyID") %>% 
  full_join(urchins, by="SurveyID") %>% 
  full_join(cti, by="SurveyID") %>% 
  full_join(sharks, by="SurveyID" ) %>% 
  left_join(hab_in[, c(1,4)], by="SurveyID")  %>% 
  distinct()

metrics_global$B20 <- if_else(is.na(metrics_global$B20), 0, as.numeric(metrics_global$B20))
metrics_global$fish.richness <- if_else(is.na(metrics_global$fish.richness), 0, as.numeric(metrics_global$fish.richness))
metrics_global$CF.richness <- if_else(is.na(metrics_global$CF.richness), 0, as.numeric(metrics_global$CF.richness))
metrics_global$invert.richness <- if_else(is.na(metrics_global$invert.richness), 0, as.numeric(metrics_global$invert.richness))
metrics_global$Urchins <- if_else(is.na(metrics_global$Urchins), 0, as.numeric(metrics_global$Urchins))
metrics_global$sharks <- if_else(is.na(metrics_global$sharks), 0, as.numeric(metrics_global$sharks))

metrics_global[, 2:3] <- NULL 

metrics_global <- metrics_global[!duplicated(metrics_global$SurveyID),] ##drop SiteCode duplications 
nrow(metrics_global) == length(unique(metrics_global$SurveyID)) ##check there are no SiteCode duplications 

#write.csv(metrics_global, "data/data-processed/metrics_global.csv")
save(metrics_global, file="data/data-processed/metrics_global.rda")



#### CREATE survey_indicator_values.csv table for input into Reef Life Explorer system  ----------------------------------------------

metrics_long <- gather(metrics_global, "indicator", "value", 8:ncol(metrics_global), na.rm=T)  ##change to 9:17 if habitatincluded in the metrics_global table above
metrics_long2 <- merge(metrics_long, metric_id, "indicator")
survey_indicator_values <- metrics_long2 %>%  
  dplyr::select(SurveyID, indicator_id, value) %>% 
  na.omit()

survey_indicator_values <- tibble::rowid_to_column(survey_indicator_values, "id") ##add an id column (necessary for data input into RHE system)
colnames(survey_indicator_values)[2] <- "survey_id"
survey_indicator_values <- survey_indicator_values[!(survey_indicator_values$indicator_id == "9"), ] ###removes LogB20 from indicator table - not using this anymore  

write.csv(survey_indicator_values, "outputs/data-outputs/survey_indicator_values.csv", row.names=FALSE) 

