###########################################################################
#------Script for prepping PADA data tables for use in temporal comparison----

###########################################################################

# --------------1. Wrangling PADA data tables-------------------------
# Load in PADA data tables
library(here)
library(dplyr)
creel_sample <- read.csv(here("PADA", "CREEL_SAMPLE_TABLE.csv"))
acid_rain <- read.csv(here("PADA","EXTERNAL_ACID_RAIN_FINALIZED_DATA.csv"))
acid_rain_charr <- read.csv(here("PADA","EXTERNAL_ARCTIC_CHAR_Y2013M01D28AM.csv"))
dfo_ext <- read.csv(here("PADA","EXTERNAL_DFO_PAPER_FILE_MATERIAL.csv"))
waldron <- read.csv(here("PADA","EXTERNAL_WALDRON_CODES_DATA_FINALIZED.csv"))
fish_ind_effort <- read.csv(here("PADA","FISHERIES_IND_EFFORT_TABLE.csv"))
fish_ind_sample <- read.csv(here("PADA","FISHERIES_IND_SAMPLE_TABLE.csv"))
species <- read.csv(here("PADA","SPECIES_LIST.csv"))
waterbodies <- read.csv(here("PADA","WATERBODIES_TABLE.csv"))

# Filter only for columns needed in each table, filter for missing data, etc.
# 1. Creel Survey Data
creel_sample <- creel_sample[,c(3:7,11,14,15)]
waterbodies <- waterbodies[,c(4:7,9,10,16:18)]
colnames(waterbodies)[4] <- "LAKE_NAME"
write.csv(waterbodies,"waterbodies2.csv",row.names=FALSE)

# Combine creel_sample and creel_effort dfs together:
creel <- merge(creel_sample, waterbodies, by=c("LAKE_NAME", "ROUTE_NAME"),
               all.x=TRUE, all.y=FALSE)
# Exclude anything without a LAKE_NAME (checked in original file, no coords either, so useless)
creel <- creel %>% filter(LAKE_NAME!="")
# Check for entries without coords
creel %>% filter(is.na(LAT_Y_COORD)) %>% distinct(LAKE_NAME)
creel %>% filter(is.na(LONG_X_COORD)) %>% distinct(LAKE_NAME)
#Two LAKE_NAMES: Lobstick and Ship harbour Little Pond; add these coords manually
creel[creel$LAKE_NAME=="Lobstick", "LAT_Y_COORD"] <- 53.86346
creel[creel$LAKE_NAME=="Lobstick", "LONG_X_COORD"] <- -64.94744
creel[creel$LAKE_NAME=="Ship harbour Little Pond", "LAT_Y_COORD"] <- 47.487688
creel[creel$LAKE_NAME=="Ship harbour Little Pond", "LONG_X_COORD"] <- -53.706088
# Merge with species data:
species <- species[,1:3]
creel_spec <- merge(creel, species, by="SPEC_CODE",all.x=TRUE, all.y=FALSE)
# Exclude anything without a species common name (n=5)
creel_spec <- creel_spec %>% filter(COMMON_NAME!="") %>% filter(!is.na(COMMON_NAME))

# 2. Acid rain survey data
acid_rain <- acid_rain[,c(4,9:11,14,15,19,20,25,26)]
acid_rain_charr <- acid_rain_charr[,c(4,9:14,16,17,20)]
# bind these two dfs together, filling empty columns with NAs
colnames(acid_rain)[8] <- "SPEC_CODE"
colnames(acid_rain_charr)[4] <- "WF_NAME"
library(plyr)
acid_rain_spec <- rbind.fill(acid_rain,acid_rain_charr)
detach(package:plyr)
# Check for entries without coords
acid_rain_spec %>% filter(is.na(LAT_DD)) # only one; seems to be a test entry; therefore exclude
acid_rain_spec %>% filter(is.na(LONG_DD))
acid_rain_spec <- acid_rain_spec %>% filter(!is.na(LAT_DD))
# Merge with species data:
acid_rain_spec <- merge(acid_rain_spec, species, by="SPEC_CODE",all.x=TRUE, all.y=FALSE)
# Everything has a species common name

# Load in required package for matching closest coordinates
# want to match to closest waterbodies entry to get approximate ecoregion data (will still use original coords in analysis)
library(hutilscpp)
acid_rain_match <- match_nrst_haversine(acid_rain_spec$LAT_DD, 
                                        acid_rain_spec$LONG_DD,
                                        waterbodies$LAT_Y_COORD,
                                        waterbodies$LONG_X_COORD,
                                        close_enough="0.25km",
                                        .verify_box = TRUE)
# Check how many entries have a match greater than 1km away to waterbodies:
acid_rain_match %>% filter(dist>1)
max(acid_rain_match$dist)
# furthest distance is 12.25km
# bind back with original hist_coords df (in same order as started):
acid_rain_match <- cbind(acid_rain_spec, acid_rain_match)
# Bind with matching entry in waterbodies df
waterbodies$pos <- as.integer(row.names(waterbodies))
acid_rain_merge <- merge(acid_rain_match, waterbodies, by="pos")


# 3. DFO external paper file material
dfo_ext <- dfo_ext[,c(3,4,10,15,16,18:20,23,26,30,33,47)]
dfo_ext %>% filter(is.na(LAT_DD)) # No lat coords missing
dfo_ext %>% filter(is.na(LONG_DD)) %>% distinct(LAKE)
# Many longs missing: Lobstick, Mobile Big Pond, Char Lake, Goose Pond, 
# Hardings Pond, Long Pond, Soldiers Pond, Southwest Pond, Three Mile Lake
dfo_ext[dfo_ext$LAKE=="Lobstick", "LAT_DD"] <- 53.86346
dfo_ext[dfo_ext$LAKE=="Lobstick", "LONG_DD"] <- -64.94744
dfo_ext[dfo_ext$LAKE=="Mobile Big Pond", "LAT_DD"] <- 47.268000001
dfo_ext[dfo_ext$LAKE=="Mobile Big Pond", "LONG_DD"] <- -53.01500001
dfo_ext[dfo_ext$LAKE=="Char Lake", "LAT_DD"] <- 58.210000001
dfo_ext[dfo_ext$LAKE=="Char Lake", "LONG_DD"] <- -63.04300001
dfo_ext[dfo_ext$LAKE=="Goose Pond", "LAT_DD"] <- 47.43293
dfo_ext[dfo_ext$LAKE=="Goose Pond", "LONG_DD"] <- -53.51857
dfo_ext[dfo_ext$LAKE=="Hardings Pond", "LAT_DD"] <- 49.56492
dfo_ext[dfo_ext$LAKE=="Hardings Pond", "LONG_DD"] <- -57.74817
dfo_ext[dfo_ext$LAKE=="Long Pond", "LAT_DD"] <- 47.3485
dfo_ext[dfo_ext$LAKE=="Long Pond", "LONG_DD"] <- -52.81226
dfo_ext[dfo_ext$LAKE=="Soldiers Pond", "LAT_DD"] <- 47.407000001
dfo_ext[dfo_ext$LAKE=="Soldiers Pond", "LONG_DD"] <- -52.994000001
dfo_ext[dfo_ext$LAKE=="Southwest Pond", "LAT_DD"] <- 47.34806
dfo_ext[dfo_ext$LAKE=="Southwest Pond", "LONG_DD"] <- -53.22593
dfo_ext[dfo_ext$LAKE=="Three Mile Lake", "LAT_DD"] <- 50.982000001
dfo_ext[dfo_ext$LAKE=="Three Mile Lake", "LONG_DD"] <- -56.871000001

# merge with species
colnames(dfo_ext)[10] <- "SPEC_CODE"
dfo_ext_spec <- merge(dfo_ext, species, by="SPEC_CODE",all.x=TRUE, all.y=FALSE)
# Everything has a species common name

# Match with nearest location in waterbodies table to get ecoregion data:
dfo_ext_match <- match_nrst_haversine(dfo_ext_spec$LAT_DD, 
                                      dfo_ext_spec$LONG_DD,
                                      waterbodies$LAT_Y_COORD,
                                      waterbodies$LONG_X_COORD,
                                      close_enough="0.25km",
                                      .verify_box = TRUE)
# Check how many entries have a match greater than 1km away to waterbodies:
dfo_ext_match %>% filter(dist>1)
max(dfo_ext_match$dist)
# furthest distance is 14.9km
# bind back with original hist_coords df (in same order as started):
dfo_ext_match <- cbind(dfo_ext_spec, dfo_ext_match)
# Bind with matching entry in waterbodies df
waterbodies$pos <- as.integer(row.names(waterbodies))
dfo_ext_merge <- merge(dfo_ext_match, waterbodies, by="pos")

# 4. External waldron codes
waldron <- waldron[,c(5,10,15:18,20,22,23,25:27)]
# One row with no LAT_DD, seems empty completely- therefore exclude
waldron <- waldron %>% filter(!is.na(LAT_DD))
waldron <- waldron %>% filter(!is.na(LONG_DD))
colnames(waldron)[7] <- "SPEC_CODE"
# merge with species
waldron_spec <- merge(waldron, species, by="SPEC_CODE",all.x=TRUE, all.y=FALSE)
# Everything has a species common name

# Match with nearest location in waterbodies table to get ecoregion data:
waldron_match <- match_nrst_haversine(waldron_spec$LAT_DD, 
                                      waldron_spec$LONG_DD,
                                      waterbodies$LAT_Y_COORD,
                                      waterbodies$LONG_X_COORD,
                                      close_enough="0.25km",
                                      .verify_box = TRUE)
# Check how many entries have a match greater than 1km away to waterbodies:
waldron_match %>% filter(dist>1)
max(waldron_match$dist)
# furthest away is 5.7km
# bind back with original hist_coords df (in same order as started):
waldron_match <- cbind(waldron_spec, waldron_match)
# Bind with matching entry in waterbodies df
waterbodies$pos <- as.integer(row.names(waterbodies))
waldron_merge <- merge(waldron_match, waterbodies, by="pos")

# 5. Fisheries Independent Samples
fish_ind_sample <- fish_ind_sample[,c(4:9,12,16,17)]
# Merge with waterbodies coord data
colnames(fish_ind_sample)[3] <- "LAKE_NAME"
colnames(fish_ind_sample)[2] <- "ROUTE_NAME"
fish_ind_spec <- merge(fish_ind_sample,waterbodies,by=c("LAKE_NAME","ROUTE_NAME"),
                       all.x=TRUE, all.y=FALSE)
# Check for entries without coords
fish_ind_spec %>% filter(is.na(LAT_Y_COORD))
fish_ind_spec %>% filter(is.na(LONG_X_COORD))
#Two LAKE_NAMES: test record (exclude), and Four Mile pond
fish_ind_spec <- fish_ind_spec %>% filter(LAKE_NAME!="Donnie Test Record")
fish_ind_spec[fish_ind_spec$LAKE_NAME=="Four Mile pond", "LAT_Y_COORD"] <- 47.33727
fish_ind_spec[fish_ind_spec$LAKE_NAME=="Four Mile pond", "LONG_X_COORD"] <- -53.10599
# Merge with species data:
fish_ind_spec <- merge(fish_ind_spec, species, by="SPEC_CODE",all.x=TRUE, all.y=FALSE)
# Everything has a species common name


# Black et al. 1986 data:
historic <- read.csv(here("R output files","historic_dist_data_Lab.csv"))
# Exclude 30 locations with no coord data:
historic <- historic %>% filter(!is.na(Lat))


# Match with nearest location in waterbodies table to get ecoregion data:
hist_match <- match_nrst_haversine(historic$Lat, 
                                   historic$Long,
                                   waterbodies$LAT_Y_COORD,
                                   waterbodies$LONG_X_COORD,
                                   close_enough="0.25km",
                                   .verify_box = TRUE)
# Check how many entries have a match greater than 1km away to waterbodies:
hist_match %>% filter(dist>1)
max(hist_match$dist)
# furthest away is 34.07km
# bind back with original hist_coords df (in same order as started):
hist_match <- cbind(historic, hist_match)
# Bind with matching entry in waterbodies df
waterbodies$pos <- as.integer(row.names(waterbodies))
hist_merge <- merge(hist_match, waterbodies, by="pos")

# Have 6 complete DFs now: creel_spec, acid_rain_merge, dfo_ext_merge, 
# waldron_merge, fish_ind_spec, hist_merge (Black et al.1986)

# Combine dfs into one:
creel_toplot <- creel_spec %>%
  select(LAKE_NAME,YEAR,LAT_Y_COORD,LONG_X_COORD,SCIENTIFIC_NAME,
         ECOREGION_NAME,SUB_REGION_NAME) %>%
  mutate(source="creel") %>%
  rename(location=LAKE_NAME,year=YEAR,lat=LAT_Y_COORD,
         long=LONG_X_COORD,spec_name=SCIENTIFIC_NAME,
         ecoreg=ECOREGION_NAME, subreg=SUB_REGION_NAME)

acidrain_toplot <- acid_rain_merge %>% 
  select(WF_NAME,YEAR,LAT_DD,LONG_DD,SCIENTIFIC_NAME,
         ECOREGION_NAME,SUB_REGION_NAME) %>%
  mutate(source="gillnet") %>%
  rename(location=WF_NAME,year=YEAR,lat=LAT_DD,
         long=LONG_DD,spec_name=SCIENTIFIC_NAME,
         ecoreg=ECOREGION_NAME,subreg=SUB_REGION_NAME)

dfo_toplot <- dfo_ext_merge %>%
  mutate(source=case_when(GEAR_TYP=="" ~ "unknown",
                          GEAR_TYP=="Bottom Gill Net" ~ "gillnet",
                          GEAR_TYP=="Eel Pot"~"pot",
                          GEAR_TYP=="Electroshocker or dip net" ~"electrofish",
                          GEAR_TYP=="Fly"~ "rod",
                          GEAR_TYP=="Fyke Nets"~ "fykenet",
                          GEAR_TYP=="Gill Net"~ "gillnet",
                          GEAR_TYP=="Gill Nets"~"gillnet",
                          GEAR_TYP=="Hoop net"~"fykenet",
                          GEAR_TYP=="Hoop net & Fyke Net"~"fykenet",
                          GEAR_TYP=="Lure & Bait"~"rod",
                          GEAR_TYP=="SEINE"~"seine",
                          GEAR_TYP=="Stream trap or Weir"~"streamtrap_weir",
                          GEAR_TYP=="Surface Gill Net"~"gillnet",
                          GEAR_TYP=="Trapnet"~"fykenet",
                          GEAR_TYP=="Angling"~"rod",
                          GEAR_TYP=="Dip Net"~"dipnet",
                          GEAR_TYP=="Electro Fished"~"electrofish",
                          GEAR_TYP=="Electroshocker or Dip net"~"electrofish",
                          GEAR_TYP=="Fyke Net"~"fykenet",
                          GEAR_TYP=="Fyke Nets & Hoop"~"fykenet",
                          GEAR_TYP=="Gill nets"~"gillnet",
                          GEAR_TYP=="Gillnet"~"gillnet",
                          GEAR_TYP=="Hoop Net"~"fykenet",
                          GEAR_TYP=="Live Trap"~"livetrap",
                          GEAR_TYP=="Not known"~"unknown",
                          GEAR_TYP=="Pots (Baited0"~"pot",
                          GEAR_TYP=="ROD"~"rod",
                          GEAR_TYP=="Seine Net"~"seine",
                          GEAR_TYP=="Stream Trap or Weir"~"streamtrap_weir",
                          GEAR_TYP=="Trap Nets & Gill Nets"~"gillnet")) %>%
  select(LAKE,LAT_DD,LONG_DD,YEAR,SCIENTIFIC_NAME,
         ECOREGION_NAME,SUB_REGION_NAME,source) %>%
  rename(location=LAKE,year=YEAR,lat=LAT_DD,
         long=LONG_DD,spec_name=SCIENTIFIC_NAME,
         ecoreg=ECOREGION_NAME,subreg=SUB_REGION_NAME)

waldron_toplot <- waldron_merge %>% 
  select(LOCATION,LAT_DD,LONG_DD,YEAR,SCIENTIFIC_NAME,
         ECOREGION_NAME,SUB_REGION_NAME) %>%
  mutate(source="unknown") %>%
  rename(location=LOCATION,year=YEAR,lat=LAT_DD,
         long=LONG_DD,spec_name=SCIENTIFIC_NAME,
         ecoreg=ECOREGION_NAME,subreg=SUB_REGION_NAME)

fish_ind_toplot <- fish_ind_spec %>% 
  mutate(source=case_when(GEAR_TYPE==""~"unknown",
                          GEAR_TYPE=="1"~"unknown",
                          GEAR_TYPE=="Angling"~"rod",
                          GEAR_TYPE=="Electrofisher"~"electrofish",
                          GEAR_TYPE=="Fyke Net"~"fykenet",
                          GEAR_TYPE=="Gill Net"~"gillnet",
                          GEAR_TYPE=="Rod"~"rod")) %>%
  select(LAKE_NAME,YEAR,LAT_Y_COORD,LONG_X_COORD,SCIENTIFIC_NAME,
         ECOREGION_NAME,SUB_REGION_NAME,source) %>%
  rename(location=LAKE_NAME,year=YEAR,lat=LAT_Y_COORD,
         long=LONG_X_COORD,spec_name=SCIENTIFIC_NAME,
         ecoreg=ECOREGION_NAME,subreg=SUB_REGION_NAME)

hist_toplot <- hist_merge %>%
  select(Location.or.Area,Lat, Long, taxon, source, 
         ECOREGION_NAME,SUB_REGION_NAME) %>%
  rename(location=Location.or.Area, lat=Lat, long=Long,
         spec_name=taxon, ecoreg=ECOREGION_NAME,
         subreg=SUB_REGION_NAME) %>%
  mutate(year=case_when(source=="Blacketal" ~ "1986",
                        source=="Perry&Keefe" ~ "2021",
                        source=="Footprints" ~ "1977",
                        source=="Michaudetal" ~ "2011")) %>%
  mutate(source=case_when(source=="Blacketal"~"unknown",
                          source=="Michaudetal"~"gillnet",
                          TRUE~ source)) %>%
  filter(source=="unknown") #exclude others for now, can plot later

# Combine all together
pada <- rbind(creel_toplot,acidrain_toplot,dfo_toplot,
              waldron_toplot,fish_ind_toplot,hist_toplot)

# Some years are coded as just '85' (for example); change to be four digits
conv_year <- function(x) {
  ifelse(x >10 & x <100, x+1900, x)
}

pada$year <- lapply(pada$year, conv_year) %>% as.numeric()
pada$year <- lapply(pada$year, conv_year) %>% as.numeric() 
# make sure it worked:
pada %>% filter(year<100)
# Note: an additional 12 samples have years coded as "1" "0" etc. or NA
# These will be excluded downstream (not useful, year data is not informative)

# Check that there isn't overlap in taxa names due to typos in database etc:
# Change Esox Lucius to Esox lucius
# Change Cottus bairdi to Cottus bairdii
# Change Catostomus commeroni to Catostomus commersonii
pada <- pada %>% 
  mutate(spec_name=case_when(spec_name=="Catostomus commersoni" ~ "Catostomus commersonii",
                             spec_name=="Cottus bairdi" ~ "Cottus bairdii",
                             spec_name=="Esox Lucius" ~ "Esox lucius",
                             spec_name=="Fundulus diaphanous" ~ "Fundulus diaphanus",
                             TRUE ~ spec_name))

# There are some pada entries with no source
pada <- pada %>%
  mutate(source=case_when(source= is.na(source) ~ "unknown",
                          TRUE ~ source))

# save full pada dataframe:
write.csv(pada,"pada_full.csv",row.names=FALSE)

# ----Only keep entries with "stream", "creek", ----
# "brook", "race", or "river" in the location name
# Also exclude "lake", "pond"
# this is to ensure best possible matching of habitat type for environmental association analysis
pada_river <- pada %>% filter(grepl("Stream", location) | 
                                grepl("stream",location) |
                                grepl("Creek",location) |
                                grepl("creek", location) |
                                grepl("Brook", location) |
                                grepl("brook", location) |
                                grepl("Race", location) |
                                grepl("race", location) |
                                grepl("River", location) |
                                grepl("river", location)) %>%
  filter(!grepl("Lake",location)) %>%
  filter(!grepl("lake",location)) %>%
  filter(!grepl("Pond",location)) %>%
  filter(!grepl("pond",location))
# This keeps 45593 entries out of 177984 (25.6%)

# Separate dataset into regions described in Black et al. 1986 paper,
# Add a column for this region and bind back together
lakemelpada_river <- pada_river %>% filter(ecoreg=="High Boreal Forest (Lake Melville)") %>%
  mutate(fish_reg="Lake Melville")
southeastpada_river <- pada_river %>% filter(ecoreg!="High Boreal Forest (Lake Melville)" &
                                               ecoreg!="Avalon Forest" &
                                               ecoreg!="Central Newfoundland Forest" &
                                               ecoreg!="Eastern Hyper-Oceanic Barrens" &
                                               ecoreg!="Long Range Barrens" &
                                               ecoreg!="Long Range Barrens (Questionable Ecoregion)" &
                                               ecoreg!="Maritime Barrens" &
                                               ecoreg!="Maritime Barrens (Questionable Ecoregion" &
                                               ecoreg!="North Shore Forest" &
                                               ecoreg!="Northern Peninsula Forest" &
                                               ecoreg!="Strait Of Belle Isle" &
                                               ecoreg!="Western Newfoundland Forest") %>%
  filter(lat<54.23333) %>%
  mutate(fish_reg="Southeast Labrador")
northfraserpada_river <- pada_river %>% filter(ecoreg!="High Boreal Forest (Lake Melville)") %>%
  filter(lat>54.23333 & lat<56.61884) %>%
  mutate(fish_reg="Northern Labrador (to Fraser River)")
northnorthpada_river <- pada_river %>% filter(lat>56.61884) %>%
  mutate(fish_reg="Northern Labrador (above Fraser River)")
nfldpada_river <- pada_river %>%
  filter(ecoreg=="Avalon Forest" |
           ecoreg=="Central Newfoundland Forest" |
           ecoreg=="Eastern Hyper-Oceanic Barrens" |
           ecoreg=="Long Range Barrens" |
           ecoreg=="Long Range Barrens (Questionable Ecoregion)" |
           ecoreg=="Maritime Barrens" |
           ecoreg=="Maritime Barrens (Questionable Ecoregion" |
           ecoreg=="North Shore Forest" |
           ecoreg=="Northern Peninsula Forest" |
           ecoreg=="Strait Of Belle Isle" |
           ecoreg=="Western Newfoundland Forest") %>%
  mutate(fish_reg="Newfoundland")
pada_river2<- rbind(lakemelpada_river,southeastpada_river,
                    northfraserpada_river,northnorthpada_river,
                    nfldpada_river)
write.csv(pada_river2,"pada_river.csv",row.names=FALSE)

