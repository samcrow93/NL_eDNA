###########################################################################
#------Script for environmental association analyses of eDNA community composition data----

###########################################################################

#-------Prepping eDNA data for PCoA------------
# load in eDNA data
# also load in working waterbodies2.csv file (from PADA dataset, generated 4)_Prepping_PADA.R)
library(here)
library(dplyr)
cert3_geog <- read.csv(here("R output files","cert3_geog.csv"))
waterbodies <- read.csv(here("R output files","waterbodies2.csv"))
# Format eDNA data to leave out non-species level
# Also drop marine species: M villosus, Ammodytes, B. saida, E. praecisus,
# L. gibbus, Myoxocephalus, G. morhua, P. americanus, T. adspersus
cert3_spec <- cert3_geog %>% filter(grepl(" ",taxon)) %>%
  filter(taxon!="Mallotus villosus" &
           taxon!="Boreogadus saida" & taxon!="Eumesogrammus praecisus" &
           taxon!="Liparis gibbus" & taxon!="Gadus morhua" &
           taxon!="Pseudopleuronectes americanus" & taxon!="Tautogolabrus adspersus")

edna <- cert3_spec %>% filter(prop_site_consensus>0) %>%
  select(RiverCode,year, Latitude, Longitude, taxon) %>%
  rename(location=RiverCode, lat=Latitude, long=Longitude,
         spec_name=taxon) %>%
  mutate(source="eDNA")

# Match to closest location in waterbodies df to get approximate ecoregion and subregion
library(hutilscpp)
edna_match <- match_nrst_haversine(edna$lat, 
                                   edna$long,
                                   waterbodies$LAT_Y_COORD,
                                   waterbodies$LONG_X_COORD,
                                   close_enough="0.25km",
                                   .verify_box = TRUE)
# Check how many entries have a match greater than 1km away to waterbodies:
edna_match %>% filter(dist>1)
max(edna_match$dist)
hist(edna_match$dist)
# There are quite a few, furthest is 46.7km away, but most are <20km away
# bind back with original hist_coords df (in same order as started):
edna_match <- cbind(edna, edna_match)
# Bind with matching entry in waterbodies df
waterbodies$pos <- as.integer(row.names(waterbodies))
edna_merge <- merge(edna_match, waterbodies, by="pos")

edna2 <- edna_merge %>%
  select(location, year,lat,long,spec_name,source,ECOREGION_NAME,SUB_REGION_NAME) %>%
  rename(ecoreg=ECOREGION_NAME,subreg=SUB_REGION_NAME)

# Add in column for regions described in Black et al. 1986 paper:
lakemel_edna <- edna2 %>% filter(ecoreg=="High Boreal Forest (Lake Melville)") %>%
  mutate(fish_reg="Lake Melville")
southeast_edna <- edna2 %>% filter(ecoreg!="High Boreal Forest (Lake Melville)" &
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
northfraser_edna <- edna2 %>% filter(ecoreg!="High Boreal Forest (Lake Melville)") %>%
  filter(lat>54.23333 & lat<56.61884) %>%
  mutate(fish_reg="Northern Labrador (to Fraser River)")
northnorth_edna <- edna2 %>% filter(lat>56.61884) %>%
  mutate(fish_reg="Northern Labrador (above Fraser River)")
nfld_edna <- edna2 %>%
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
edna3<- rbind(lakemel_edna,southeast_edna,
              northfraser_edna,northnorth_edna,
              nfld_edna)
write.csv(edna3,"edna3.csv", row.names=FALSE) # save for use with PADA comparisons
edna4 <- edna3[,-10] #get rid of pos column now
edna4 <- edna4 %>% select(location,lat,long,
                          spec_name,source,fish_reg) %>%
  rename(RiverCode=location)

# Convert edna4 to wide:
library(tidyr)
edna_wide2 <- edna4 %>%
  group_by(lat, long,RiverCode,spec_name,source,fish_reg) %>%
  summarise(n_obs=n()) %>%
  spread(spec_name,n_obs) %>%
  as.data.frame() %>%
  replace(is.na(.),0)
# Convert observation counts to just 1 (pres/abs format)
edna_wide2 <- as.data.frame(cbind(edna_wide2[,1:5],(+(edna_wide2[-(1:5)] > 0))))

# Check for sites with no data now:
rowSums(edna_wide2[,6:31]) #all rows (sites) have data
# Arrange by RiverCode
edna_wide2 <- edna_wide2 %>% arrange(factor(RiverCode))
write.csv(edna_wide2,"edna_wide2.csv",row.names=FALSE)

#--------PCoA of all eDNA sites-----------
library(vegan)
pcoa_edna <- capscale(edna_wide2[,c(6:31)] ~ 1, dist="jaccard")
smry_pcoa_edna <- summary(pcoa_edna)
pcoa_edna_sites <- cbind(edna_wide2[,c(1:5)], 
                         data.frame(smry_pcoa_edna$sites[,1:2])) # loadings for RDA1 and 2 for each site
pcoa_edna_species <- data.frame(smry_pcoa_edna$species[,1:2]) # loadings for RDA1 and 2 for each taxon

#Biplot 
#reorder fish regs levels
colnames(pcoa_edna_sites)[5] <- "Region"
pcoa_edna_sites$Region <- factor(pcoa_edna_sites$Region, 
                                 levels=c("Newfoundland",
                                          "Southeast Labrador",
                                          "Lake Melville",
                                          "Northern Labrador (to Fraser River)",
                                          "Northern Labrador (above Fraser River)"))

pcoa_edna_plot <- ggplot(pcoa_edna_sites,
                         aes(x=MDS1,y=MDS2,color=Region))+
  geom_point(size=3)+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted") +
  coord_fixed()+
  #stat_ellipse(aes(x=MDS1,y=MDS2, group=source),type="norm")+
  scale_color_manual(values=c("Newfoundland"="maroon",
                              "Southeast Labrador"="darkorange",
                              "Lake Melville"="red",
                              "Northern Labrador (to Fraser River)" ="blue",
                              "Northern Labrador (above Fraser River)"="darkgreen"))+
  #stat_ellipse(aes(x = MDS1,y=MDS2,fill=fish_reg,group=fish_reg),
  #  geom="polygon",level=0.8,alpha=0.1)+
  #scale_fill_manual(values=c("Lake Melville"="red",
  #                          "Newfoundland"="maroon",
  #                         "Northern Labrador (above Fraser River)"="darkgreen",
  #                        "Northern Labrador (to Fraser River)" ="blue",
  #                       "Southeast Labrador"="darkorange"))+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(axis.text.x = element_text(vjust = 0.5, hjust=1,size=18,colour="black"))+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  theme(axis.title=element_text(size=18))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=18))+ #change legend text font size
  xlab("MDS1") + ylab("MDS2")
# facet_wrap(~fish_reg)+
# ggtitle("PCoA eDNA community composition data")
pcoa_edna_plot

# Find species that load the strongest on the first two PCoA axes (top 10%):
pcoa_edna_speciesb <- pcoa_edna_species %>% arrange(MDS1)
obs <- nrow(pcoa_edna_speciesb)
pcoa_edna_speciesb %>% filter(row_number() < obs * 0.1)
# top 10% most pos PCoA1 were S. alpinus and C. catostomus
pcoa_edna_speciesc <- pcoa_edna_species %>% arrange(desc(MDS1))
obsc <- nrow(pcoa_edna_speciesc)
pcoa_edna_speciesc %>% filter(row_number() < obs * 0.1) 
# top 10% most neg PCoA1 were A. rostrata and S. salar

pcoa_edna_speciesd <- pcoa_edna_species %>% arrange(MDS2)
obsd <- nrow(pcoa_edna_speciesd)
pcoa_edna_speciesd %>% filter(row_number() < obs * 0.1)
# top 10% most pos PCoA1 were G. aculeatus and S. salar
pcoa_edna_speciesd <- pcoa_edna_species %>% arrange(desc(MDS2))
obsc <- nrow(pcoa_edna_speciesd)
pcoa_edna_speciesd %>% filter(row_number() < obs * 0.1) 
# top 10% most neg PCoA1 were S. alpinus and A. rostrata

# --------------Generating Moran's Eigenvector Maps (MEMs) for use in RDAs----------
# Calculate in-water depths using marmap
library(data.table)  
library(tidyverse)
library(marmap)
library(parallel)
coords <- edna_wide2 %>% select(long, lat)

bathydata2 <- getNOAA.bathy(-68,-50,46,62, res=2,keep=T)
trans1 <- marmap::trans.mat(bathydata2)

# Compute distance from each eDNA site to nearest point at the coast
# isobath=0 is nearest point on coast, so use isobath= -10
mar_dist <- dist2isobath(bathydata2, coords$long, coords$lat, isobath = -10)
# Check depths of new coordinates:
onland <- get.depth(bathydata2, mar_dist[,4:5], locator = FALSE) %>%
  filter(depth>=0) # 28 points still have a positive depth:

notland <- get.depth(bathydata2, mar_dist[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start <- mar_dist %>% filter(end.lat %in% notland$lat)

# run dist2isobath again for these 28 points:
land <- mar_dist %>% filter(end.lat %in% onland$lat) # filter for starting coords that were still not on land
mar_dist2 <- dist2isobath(bathydata2, land$start.lon, land$start.lat, isobath = -20)
onland2 <- get.depth(bathydata2, mar_dist2[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # 14 of them are still on land
# add sites that are not not on land to notland df:
notland2 <- get.depth(bathydata2, mar_dist2[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start2 <- mar_dist2 %>% filter(end.lat %in% notland2$lat)
notland_start <- rbind(notland_start,notland_start2)

# Bump up isobath again:
land2 <- mar_dist2 %>% filter(end.lat %in% onland2$lat) # filter for starting coords that were still not on land
mar_dist3 <- dist2isobath(bathydata2, land2$start.lon, land2$start.lat, isobath = -30)
onland3 <- get.depth(bathydata2, mar_dist3[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # 7 still on land, bump up
# add sites that are not not on land to notland df:
notland3 <- get.depth(bathydata2, mar_dist3[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start3 <- mar_dist3 %>% filter(end.lat %in% notland3$lat)
notland_start <- rbind(notland_start,notland_start3)

# Bump up isobath again:
land3 <- mar_dist3 %>% filter(end.lat %in% onland3$lat) # filter for starting coords that were still not on land
mar_dist4 <- dist2isobath(bathydata2, land3$start.lon, land3$start.lat, isobath = -50)
onland4 <- get.depth(bathydata2, mar_dist4[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # 3 still on land
# add sites that are not not on land to notland df:
notland4 <- get.depth(bathydata2, mar_dist4[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start4 <- mar_dist4 %>% filter(end.lat %in% notland4$lat)
notland_start <- rbind(notland_start,notland_start4)

# Bump up isobath again:
land4 <- mar_dist4 %>% filter(end.lat %in% onland4$lat) # filter for starting coords that were still not on land
mar_dist5 <- dist2isobath(bathydata2, land4$start.lon, land4$start.lat, isobath = -70)
land5 <- get.depth(bathydata2, mar_dist5[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # None are on land anymore
# add sites that are not not on land to notland df:
notland5 <- get.depth(bathydata2, mar_dist5[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start5 <- mar_dist5 %>% filter(end.lat %in% notland5$lat)
notland_start <- rbind(notland_start,notland_start5)
# notland_start contains the starting site coordinates, as well as 
# the coordinates that are the closest in-ocean coordinates

notland_start <- notland_start %>% arrange(start.lat) #arrange by starting latitude to match up with community data later
marine_sites <- notland_start[,4:5] 
lcp <- lc.dist(trans1, marine_sites , res = "dist", meters=T)
lcp # distance matrix

# Calculate MEMs by weighting Euclidean distances by inverse of marine distance
library(dplyr)
coords <- coords %>% arrange(lat)
nb = spdep::graph2nb(spdep::gabrielneigh(coords), sym = TRUE)  # Gabriel graph: neighbor definition (a list of neighbors for each site/individual)
spdep::plot.nb(nb, coords, pch = 21)
disttri = spdep::nbdists(nb, coords,longlat=TRUE)  # Add longlat=T for lat/long coordinates
# modify the weights based on the actual geographic (marine) distances between neighbors
im = 1/lcp  # Use inverse distance weights, (the further away, the lower the value)
im[!is.finite(im)]<- 0
im2 <- as.matrix(im)

# get a revised spatial weights matrix
fdist = lapply(1:length(nb), function(x){c(sapply(nb[x], 
                                                  function(z){ im2[x,z]}))}) # end lapply
listW = spdep::nb2listw(nb, style="W", glist = fdist) #glist adds revised weights
summary(sapply(listW$weights, sum))
mems = adespatial::scores.listw(listW, MEM.autocor = "positive") 
summary(mems)
# Generates 83 positive MEMs

# Next, check for associations between MEMs and climate data
# Goal to pick uncorrelated MEMs that are best predictors
# Load in climate average data and MEM data:
all_climate2 <- read.csv(here("R output files","climate_avgs.csv"))

# merge community data with climate data:
comm_env_edna <- merge(edna_wide2, all_climate2,by=c("RiverCode","source"),
                       all.x=TRUE,all.y=FALSE)
comm_env_edna <- comm_env_edna %>% arrange(lat)
# merge also with dbmems:
comm_env_mem_edna <- cbind(comm_env_edna,mems)

# Check for correlations among explanatory variables (MEMs and climate)
cors <- data.frame(cor(comm_env_mem_edna[,32:121]))
#First check if any MEMs are strongly correlated (>0.75) with climate vars:
cors_clim <- cors[1:7,8:90]
cors_clim %>% filter_all(any_vars(abs(.) > abs(0.75)))
# There are no MEMs correlated with climate variables above 0.75

# Run RDA wit just MEMs and use ordistep to pick a subset for use in RDA
full_rda_mems <- comm_env_mem_edna[,39:121]
library(vegan)
dbrda_edna_mems <- capscale(comm_env_mem_edna[,c(6:31)] ~ .,
                            full_rda_mems,
                            dist="jaccard")
summary(dbrda_edna_mems)
# this rda explains 68.14% of the variation in the data
#https://www.nature.com/articles/s41437-020-0352-6

# Compare with null rda, which is just a pcoa
dbrda_edna_null <- capscale(comm_env_mem_edna[,c(6:31)] ~ 1,
                            full_rda_mems,
                            dist="jaccard")
selection <- ordiR2step(dbrda_edna_null, scope = formula(dbrda_edna_mems), 
                        direction="both")
selection$anova  # this gives significant MEMs, 28 in total
# save selected MEMs 
marine_edna_mems <- data.frame(cbind(full_rda_mems$MEM3,full_rda_mems$MEM4,full_rda_mems$MEM18,
                                     full_rda_mems$MEM2,full_rda_mems$MEM9,full_rda_mems$MEM10,
                                     full_rda_mems$MEM14,full_rda_mems$MEM8,full_rda_mems$MEM17,
                                     full_rda_mems$MEM6,full_rda_mems$MEM22,full_rda_mems$MEM11,
                                     full_rda_mems$MEM77,full_rda_mems$MEM79,full_rda_mems$MEM21,
                                     full_rda_mems$MEM23,full_rda_mems$MEM16,full_rda_mems$MEM56,
                                     full_rda_mems$MEM29,full_rda_mems$MEM1,full_rda_mems$MEM31,
                                     full_rda_mems$MEM49,full_rda_mems$MEM7, full_rda_mems$MEM58,
                                     full_rda_mems$MEM30,full_rda_mems$MEM78,full_rda_mems$MEM66,
                                     full_rda_mems$MEM41))
colnames(marine_edna_mems) <- c("MEM3","MEM4","MEM18",
                                "MEM2","MEM9","MEM10",
                                "MEM14","MEM8","MEM17",
                                "MEM6","MEM22","MEM11",
                                "MEM77","MEM79","MEM21",
                                "MEM23","MEM16","MEM56",
                                "MEM29","MEM1","MEM31",
                                "MEM49","MEM7", "MEM58",
                                "MEM30","MEM78","MEM66",
                                "MEM41")
write.csv(marine_edna_mems,"marine_edna_mems.csv",row.names=FALSE) # MEMs are still arranged by latitude

# --------------dbRDA with climate variables and MEMs corrected for marine distance---------
library(here)
library(dplyr)
library(vegan)
marine_mems <- read.csv(here("R output files","marine_edna_mems.csv"))
all_climate2 <- read.csv(here("R output files","climate_avgs.csv"))
edna_wide2 <- read.csv(here("R output files","edna_wide2.csv"))

# merge community data with climate data:
comm_env_edna <- merge(edna_wide2, all_climate2,by=c("RiverCode","source"),
                       all.x=TRUE,all.y=FALSE)
comm_env_edna <- comm_env_edna %>% arrange(lat)
# merge also with dbmems:
comm_env_mem_marine <- cbind(comm_env_edna,marine_mems)

# For variance partitioning, need to run 3 models:
# 1. full model
dbrda_edna_marinemem_full <- capscale(comm_env_mem_marine[,c(6:31)] ~ mean_max_air_temp +
                                        mean_rel_humid + mean_wind_speed_2m +
                                        mean_solar_rad + mean_atm_press +
                                        mean_snow_precip + mean_snow_accum +
                                        MEM3 + MEM4 +MEM18 +
                                        MEM2 + MEM9 + MEM10 +
                                        MEM14 +MEM8 + MEM17 +
                                        MEM6 + MEM22 +MEM11 +
                                        MEM77 + MEM79 + MEM21 +
                                        MEM23 + MEM16 +MEM56 +
                                        MEM29 + MEM1 + MEM31 +
                                        MEM49 + MEM7 + MEM58 +
                                        MEM30 + MEM78 + MEM66 +
                                        MEM41, 
                                      comm_env_mem_marine,
                                      dist="jaccard")
summary(dbrda_edna_marinemem_full)
anova.cca(dbrda_edna_marinemem_full)
anova.cca(dbrda_edna_marinemem_full,by="terms",permu=9999)
#amount of variation explained by the full model= 56.21% (residual= 43.79%)

# 2. pure environment (conditioned with MEMs)
dbrda_edna_marinemem <- capscale(comm_env_mem_marine[,c(6:31)] ~ mean_max_air_temp +
                                   mean_rel_humid + mean_wind_speed_2m +
                                   mean_solar_rad + mean_atm_press +
                                   mean_snow_precip + mean_snow_accum +
                                   Condition(MEM3 + MEM4 + MEM18 +
                                               MEM2 + MEM9 + MEM10 +
                                               MEM14 + MEM8 + MEM17 +
                                               MEM6 + MEM22 + MEM11 +
                                               MEM77 + MEM79 + MEM21 +
                                               MEM23 + MEM16 + MEM56 +
                                               MEM29 + MEM1 + MEM31 +
                                               MEM49 + MEM7 + MEM58 +
                                               MEM30 + MEM78 + MEM66 +
                                               MEM41), 
                                 comm_env_mem_marine,
                                 dist="jaccard")

smry_edna_rda_marinemem <- summary(dbrda_edna_marinemem)
smry_edna_rda_marinemem
# conditioned vars (28 significant MEMs) explain 48.68%
# constrained vars (climate) explain 7.54%
anova.cca(dbrda_edna_marinemem) # model is significant (p=0.001)
anova.cca(dbrda_edna_marinemem,by="terms", permu=9999)
anova.cca(dbrda_edna_marinemem,by="axis", permu=9999)
# all climate vars are significant at p<=0.001 
# first 3 axes are significant at p<=0.05

# 3. An RDA without MEMs as conditioning variables 
# This one isn't used in the variance partitioning (not necessary)
# BUT: useful for visualizing effects
# Represents total enviro influence (pure AND spatially-structured)
# Plot this one
dbrda_edna_marinemem2 <- capscale(comm_env_mem_marine[,c(6:31)] ~ mean_max_air_temp +
                                    mean_rel_humid + mean_wind_speed_2m +
                                    mean_solar_rad + mean_atm_press +
                                    mean_snow_precip + mean_snow_accum, 
                                  comm_env_mem_marine,
                                  dist="jaccard")

smry_edna_rda_marinemem2 <- summary(dbrda_edna_marinemem2)
smry_edna_rda_marinemem2
# climate variables explain 35.44% of the variation (same as model without other MEMs)
anova.cca(dbrda_edna_marinemem2)

# 4. Pure spatial (conditioned on environmental variables)
dbrda_edna_marinemem_purespace <- capscale(comm_env_mem_marine[,c(6:31)] ~ MEM3 + MEM4 +MEM18 +
                                             MEM2 + MEM9 + MEM10 +
                                             MEM14 +MEM8 + MEM17 +
                                             MEM6 + MEM22 +MEM11 +
                                             MEM77 + MEM79 + MEM21 +
                                             MEM23 + MEM16 +MEM56 +
                                             MEM29 + MEM1 + MEM31 +
                                             MEM49 + MEM7 + MEM58 +
                                             MEM30 + MEM78 + MEM66 +
                                             MEM41 + Condition(mean_max_air_temp +
                                                                 mean_rel_humid + mean_wind_speed_2m +
                                                                 mean_solar_rad + mean_atm_press +
                                                                 mean_snow_precip + mean_snow_accum),
                                           comm_env_mem_marine,
                                           dist="jaccard")
summary(dbrda_edna_marinemem_purespace)
anova.cca(dbrda_edna_marinemem_purespace)
# amount explained by MEMs (space) is 20.78%
# amount explained by climate (conditioned) is 35.44%
# total is 56.22%

# To summarize variance partitioning:
# Total amount of variance explained (enviro and spatial factors together): 56.21%
# Residual variance (unexplained): 43.79%
# Pure enviro: 7.54
# Pure spatial: 20.78%
# Confounded enviro/spatial: 56.21-7.54-20.78= 27.89%
# total enviro (including spatially-structured enviro influences) is 35.43%

#save site loadings for first two RDA axes (first RDA including MEMs):
site_loadings_rda_marinemem <- as.data.frame(smry_edna_rda_marinemem$sites[,1:2])
vars_loadings_rda_marinemem <- as.data.frame(smry_edna_rda_marinemem$biplot[,1:2])
rownames(vars_loadings_rda_marinemem) <- c("Max Air Temp",
                                           "Relative Humidity",
                                           "Wind Speed 2m",
                                           "Solar Radiation",
                                           "Atmospheric Pressure",
                                           "Snow Precip",
                                           "Snow Accum")
#bind back with comm_env_pc
climate_meta_rda_marinemem <- cbind(comm_env_mem_marine[,1:5],
                                    site_loadings_rda_marinemem)

climate_meta_rda_marinemem <- climate_meta_rda_marinemem %>% 
  arrange(source, RiverCode)

# Find species that load the strongest on the first two PCoA axes (top 10%):
vars_loadings_rda_marinemem %>% arrange(CAP1)
vars_loadings_rda_marinemem %>% arrange(CAP2)


# Biplot
colnames(climate_meta_rda_marinemem)[5] <- "Region"
ggplot(climate_meta_rda_marinemem,aes(x=CAP1,y=CAP2, col=Region))+
  geom_point(size=4,alpha=0.7)+
  scale_color_manual(values=c("Newfoundland"="maroon",
                              "Southeast Labrador"="darkorange",
                              "Lake Melville"="red",
                              "Northern Labrador (to Fraser River)" ="blue",
                              "Northern Labrador (above Fraser River)"="darkgreen"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted") +
  #coord_fixed() +
  #geom_segment(data=segments, aes(x=PC1_hist,xend=PC1_edna, 
  #    y=PC2_hist, yend=PC2_edna),
  # color="blue",arrow=arrow(length=unit(0.3,"cm")))+
  #stat_ellipse(aes(x = CAP1,y=CAP2,fill=fish_reg,group=fish_reg),
  #  geom="polygon",level=0.8, alpha=0.2)+
  # scale_fill_manual(values=c("Lake Melville"="red",
  #                           "Newfoundland"="maroon",
  #                          "Northern Labrador (above Fraser River)"="darkgreen",
  #                         "Northern Labrador (to Fraser River)" ="blue",
  #                        "Southeast Labrador"="darkorange"))+
ggrepel::geom_label_repel(data=vars_loadings_rda_marinemem,aes(x=CAP1,y=CAP2,label=rownames(vars_loadings_rda_marinemem)),
                          #  hjust=0.5*(1-sign(CAP1)),
                          #  vjust=0.5*(1-sign(CAP2))),
                          color="black",size=3.5,alpha=0.8, position=position_dodge(0))+
  geom_segment(data=vars_loadings_rda_marinemem, aes(x=0,xend=CAP1,y=0,yend=CAP2),
               color="black",size=1)+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(axis.text.x = element_text(size=12,colour="black"))+
  theme(axis.text.y = element_text(size=12,colour="black"))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14))+ #change legend t
  xlab("CAP1") + ylab("CAP2")
#facet_wrap(~fish_reg)+
#ggtitle("db-RDA eDNA community composition and climate with MEMs (corrected for marine distances)")


#save site loadings for first two RDA axes (second RDA without MEMs):
site_loadings_rda_marinemem2 <- as.data.frame(smry_edna_rda_marinemem2$sites[,1:2])
vars_loadings_rda_marinemem2 <- as.data.frame(smry_edna_rda_marinemem2$biplot[,1:2])
rownames(vars_loadings_rda_marinemem2) <- c("Max Air Temp",
                                            "Relative Humidity",
                                            "Wind Speed 2m",
                                            "Solar Radiation",
                                            "Atmospheric Pressure",
                                            "Snow Precip",
                                            "Snow Accum")
#bind back with comm_env_pc
climate_meta_rda_marinemem2 <- cbind(comm_env_mem_marine[,1:5],
                                     site_loadings_rda_marinemem2)

climate_meta_rda_marinemem2 <- climate_meta_rda_marinemem2 %>% 
  arrange(source, RiverCode)


# Biplot
colnames(climate_meta_rda_marinemem2)[5] <- "Region"
ggplot(climate_meta_rda_marinemem2,aes(x=CAP1,y=CAP2, col=Region))+
  geom_point(size=4, alpha=0.7)+
  scale_color_manual(values=c("Newfoundland"="maroon",
                              "Southeast Labrador"="darkorange",
                              "Lake Melville"="red",
                              "Northern Labrador (to Fraser River)" ="blue",
                              "Northern Labrador (above Fraser River)"="darkgreen"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted") +
  #coord_fixed() +
  #geom_segment(data=segments, aes(x=PC1_hist,xend=PC1_edna, 
  #    y=PC2_hist, yend=PC2_edna),
  # color="blue",arrow=arrow(length=unit(0.3,"cm")))+
  #stat_ellipse(aes(x = CAP1,y=CAP2,fill=fish_reg,group=fish_reg),
  #  geom="polygon",level=0.8, alpha=0.2)+
  # scale_fill_manual(values=c("Lake Melville"="red",
  #                           "Newfoundland"="maroon",
  #                          "Northern Labrador (above Fraser River)"="darkgreen",
  #                         "Northern Labrador (to Fraser River)" ="blue",
  #                        "Southeast Labrador"="darkorange"))+
ggrepel::geom_label_repel(data=vars_loadings_rda_marinemem2,aes(x=CAP1,y=CAP2,label=rownames(vars_loadings_rda_marinemem2)),
                          # hjust=0.5*(1-sign(CAP1)),
                          # vjust=0.5*(1-sign(CAP2))),
                          color="black",size=3.5,alpha=0.8, position=position_dodge(0))+
  geom_segment(data=vars_loadings_rda_marinemem2, aes(x=0,xend=CAP1,y=0,yend=CAP2),
               color="black",size=1)+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(axis.text.x = element_text(size=12,colour="black"))+
  theme(axis.text.y = element_text(size=12,colour="black"))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14))+ #change legend te
  xlab("CAP1") + ylab("CAP2")+
  xlim(-1.5,2)
#facet_wrap(~fish_reg)+
# ggtitle("db-RDA eDNA community composition and climate")


#-----dbRDAs incorporating sampling site metadata (i.e. local habitat factors)------
met2019 <- read.csv(here("Working metadata","2019_site_metadata_reduced.csv"))
met2020 <- read.csv(here("Working metadata","2020_site_metadata_reduced.csv"))
met2021 <- read.csv(here("Working metadata","2021_site_metadata_reduced.csv"))
# Select common variables among the 3 years and bind together
met2019 <- met2019[,c(1:6,9:11,13,15,16)]
colnames(met2019)[4] <- "Center.Latitude"
colnames(met2019)[5] <- "Center.Longitude"
met2019$year <- "2019"
met2020 <- met2020[,c(1:6,9:14)]
met2020$year <- "2020"
met2021 <- met2021[,c(1:6,9,10,12,13,16,17)]
met2021$year <- "2021"
met_all <- rbind(met2019,met2020,met2021)

# There are sites that are missing some of the metadata variables
# Need to exclude because RDAs won't work with missing data
met_complete <- met_all[complete.cases(met_all),]
# this drops 33 site x year combos (out of 190) - 157 now

# Check which variables are correlated with each other:
library(tidyverse)

# load necessary functions
remove_zero_variance_columns <- function(dat) {
  out = lapply(dat, function(x) length(unique(x)))
  z_var = names(which(!out > 1))
  return(dat[!names(dat) %in% z_var])
}

#' Take a data frame with ONLY a set of predictor variables
#' Create a new data frame with only a set of variables for which the 
#' correlation coefficients do not exceede the specified threshold (default is 0.7).
get_uncorrelated_predictors = function(meta_vals, 
                                       threshold = 0.75, drop_zero_variance = TRUE){
  if(drop_zero_variance == TRUE){
    meta_vals = remove_zero_variance_columns(meta_vals)}
  tmp = cor(meta_vals)
  tmp[upper.tri(tmp)] = 0
  diag(tmp) = 0
  meta_uncorr = meta_vals[, !apply(tmp, 2, function(x) any(abs(x) > 0.75, 
                                                           na.rm = TRUE))]
  return(meta_uncorr)
}

meta_uncorr_all <- get_uncorrelated_predictors(met_complete[,6:12])
# Keeps all of them

#Bind back together:
all_met2 <- cbind(met_complete[,c(1:5)],meta_uncorr_all)
colnames(met_complete)[4] <- "lat"
colnames(met_complete)[5] <- "long"

# Check for correlations with climate data and MEMs
# First, load in all relevant dfs and combine:
marine_mems <- read.csv(here("R output files","marine_edna_mems.csv"))
all_climate2 <- read.csv(here("R output files","climate_avgs.csv"))
edna_wide2 <- read.csv(here("R output files","edna_wide2.csv"))

# merge community data with climate data:
comm_env_edna <- merge(edna_wide2, all_climate2,by=c("RiverCode","source"),
                       all.x=TRUE,all.y=FALSE)
comm_env_edna <- comm_env_edna %>% arrange(lat)
# merge also with dbmems:
comm_env_mem_marine <- cbind(comm_env_edna,marine_mems)

# subset full df based on sites that have metadata:
subset_comm <- merge(met_complete,comm_env_mem_marine,by=c("lat","long","RiverCode"),
                     all.x=TRUE,all.y=FALSE)

# Now, finally check for correlations amongst all possible explanatory vars for these 157 locations:
subset_uncorr_all <- get_uncorrelated_predictors(subset_comm[,c(6:12,42:76)])
# Doesn't drop anything; keep subset_comm df
# However, there are a couple sites for which there is multiple metadata entries (e.g. slightly diff RiverCodes such as HUR2 vs. HUR)
# Can exclude, since amalgamating community data for each site x year combo
subset_comm <- subset_comm[complete.cases(subset_comm),] # drops 5 entries only

# Variance partitioning with climate, local env, & MEMs
# Same as with full dataset, need to run 4 RDAs:
# 1. full model:
library(vegan)
subset_RDA_full <- capscale(subset_comm[,c(16:41)] ~ mean_max_air_temp +
                              mean_rel_humid + mean_wind_speed_2m +
                              mean_solar_rad + mean_atm_press +
                              mean_snow_precip + mean_snow_accum +
                              WaterTemp + DOmg.L + C.us..cm + pH.mV +
                              ORP.mV + Mean.Depth + Mean.Flow + 
                              MEM3 + MEM4 +MEM18 +
                              MEM2 + MEM9 + MEM10 +
                              MEM14 +MEM8 + MEM17 +
                              MEM6 + MEM22 +MEM11 +
                              MEM77 + MEM79 + MEM21 +
                              MEM23 + MEM16 +MEM56 +
                              MEM29 + MEM1 + MEM31 +
                              MEM49 + MEM7 + MEM58 +
                              MEM30 + MEM78 + MEM66 +
                              MEM41, 
                            subset_comm,
                            dist="jaccard")
summary(subset_RDA_full)
anova.cca(subset_RDA_full)
# Full RDA explains 62.97% of the variation (inertia explained=24.73; total inertia=39.28)
# Unexplained = 37.03% (inertia=14.54)
# RDA is significant (F= 4.413, p=0.001)

# 2. pure climate (conditioned with MEMs and local env)
subset_rda_clim <- capscale(subset_comm[,c(16:41)] ~ mean_max_air_temp +
                              mean_rel_humid + mean_wind_speed_2m +
                              mean_solar_rad + mean_atm_press +
                              mean_snow_precip + mean_snow_accum +
                              Condition(MEM3 + MEM4 + MEM18 +
                                          MEM2 + MEM9 + MEM10 +
                                          MEM14 + MEM8 + MEM17 +
                                          MEM6 + MEM22 + MEM11 +
                                          MEM77 + MEM79 + MEM21 +
                                          MEM23 + MEM16 + MEM56 +
                                          MEM29 + MEM1 + MEM31 +
                                          MEM49 + MEM7 + MEM58 +
                                          MEM30 + MEM78 + MEM66 +
                                          MEM41 +
                                          WaterTemp + DOmg.L + C.us..cm + pH.mV +
                                          ORP.mV + Mean.Depth + Mean.Flow ), 
                            subset_comm,
                            dist="jaccard")
summary(subset_rda_clim)
# Pure climate explains 6.84% of the variation (inertia=22.042; total inertia=39.276)
# Conditioning vars explain 56.122% (inertia=22.042)
# unexplained = 37.032% (inertia= 14.544)
anova.cca(subset_rda_clim)
anova.cca(subset_rda_clim,by="terms", permu=9999)
# RDA is significant (F=2.879, p=0.001)
# all climate terms are significant at p=0.005 EXCEPT wind speed and snow precip

# 3. Pure Space (conditioned on local env and climate):
subset_rda_space <- capscale(subset_comm[,c(16:41)] ~ MEM3 + MEM4 + MEM18 +
                               MEM2 + MEM9 + MEM10 +
                               MEM14 + MEM8 + MEM17 +
                               MEM6 + MEM22 + MEM11 +
                               MEM77 + MEM79 + MEM21 +
                               MEM23 + MEM16 + MEM56 +
                               MEM29 + MEM1 + MEM31 +
                               MEM49 + MEM7 + MEM58 +
                               MEM30 + MEM78 + MEM66 +
                               MEM41 +
                               Condition(mean_max_air_temp +
                                           mean_rel_humid + mean_wind_speed_2m +
                                           mean_solar_rad + mean_atm_press +
                                           mean_snow_precip + mean_snow_accum +
                                           WaterTemp + DOmg.L + C.us..cm + pH.mV +
                                           ORP.mV + Mean.Depth + Mean.Flow), 
                             subset_comm,
                             dist="jaccard")
summary(subset_rda_space)
# Pure space explains 19.65% of the variation (inertia=7.72; total inertia=39.28)
# Conditioned vars (climate and local env) explain 43.31% (17.01)
# unexplained var= 37.03% (inertia=14.54)
anova.cca(subset_rda_space)
anova.cca(subset_rda_space,by="terms", permu=9999)
# RDA is significant (F=2.0662, p=0.001)
# some MEMs are significant, some aren't

# 4. Pure local env:
subset_rda_env <- capscale(subset_comm[,c(16:41)] ~ WaterTemp + DOmg.L + C.us..cm + 
                             pH.mV +
                             ORP.mV + Mean.Depth + Mean.Flow +
                             Condition(MEM3 + MEM4 + MEM18 +
                                         MEM2 + MEM9 + MEM10 +
                                         MEM14 + MEM8 + MEM17 +
                                         MEM6 + MEM22 + MEM11 +
                                         MEM77 + MEM79 + MEM21 +
                                         MEM23 + MEM16 + MEM56 +
                                         MEM29 + MEM1 + MEM31 +
                                         MEM49 + MEM7 + MEM58 +
                                         MEM30 + MEM78 + MEM66 +
                                         MEM41 + mean_max_air_temp +
                                         mean_rel_humid + mean_wind_speed_2m +
                                         mean_solar_rad + mean_atm_press +
                                         mean_snow_precip + mean_snow_accum), 
                           subset_comm,
                           dist="jaccard")
summary(subset_rda_env)
# pure local env. explains 5.179% of the variation (inertia=2.034; total inertia=39.276)
# conditioned vars (climate and MEMs) explain 57.789% (inertia=22.697)
# unexplained = 37.032% (inertia=14.544)

anova.cca(subset_rda_env)
anova.cca(subset_rda_env,by="terms", permu=9999)
# RDA is significant (F=2.1779, p=0.05)
# all vars are significant at p=0.001 EXCEPT ORP.mV, mean.Flow, and mean.Depth


#----Plot diversity levels by site------
# First, load in map data:
library(data.table)  
library(tidyverse)
library(marmap)
library(parallel)

bathydata2 <- getNOAA.bathy(-68,-50,46,62, res=2,keep=T)
plot(bathydata2)
map2=autoplot(bathydata2, geom=c("r", "c"), 
              colour="grey90", size=0.1) +
  #scale_fill_etopo(guide = FALSE) +
  scale_fill_gradient(low="gray90",high="gray100",guide=FALSE)+
  geom_contour(aes(z=z), 
               breaks=c(-100, -200, -500, -1000, -2000, -4000), 
               colour="grey90", size=0.01) +
  xlab("Degrees Longitude") +
  ylab("Degrees Latitude") +
  theme(legend.position = "none")

# Create dataframe that summarizes taxon counts (species only) per site per year
tax_count_depth <- edna3 %>% 
  filter(grepl(" ", spec_name)) %>% 
  group_by(location,lat,long,fish_reg,year) %>%
  summarise(tax_count_all = n_distinct(spec_name)) %>% data.frame()

# Map of species richness (taxon counts)
colnames(tax_count_depth)[6] <- "Species Richness"
divmap <- map2 + geom_point(data=tax_count_depth,
                            aes(x=long,y=lat,
                                col=`Species Richness`),
                            size=3, inherit.aes=F)+
  scale_color_viridis_b(begin=1, end=0, option="B")+
  theme(legend.position="right") +
  theme(axis.title=element_text(size=18))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=18))+ #change leg
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))
#ggtitle("Total number of taxa (species-level) detected per site per year")
divmap


# ------------Kmeans clustering of the community data----------
# create jaccard distance matrix of data:
library(vegan)
dist_edna <- vegdist(comm_env_mem_edna[,c(6:31)], method="jaccard")
# clustering, testing 1-10 clusters as optimal number:
edna_clust <- as.dendrogram(hclust(dist_edna, method="ward.D2"))

kmeans_edna <- cascadeKM(dist_edna, 1, 10, iter = 5000)
plot(kmeans_edna, sortg = TRUE, grpmts.plot = TRUE)
grps <- as_tibble(kmeans_edna$partition[,2]) 
colnames(grps)[1] <- "Cluster"
# Merge back with coordinate data
cluster_coords <- cbind(comm_env_mem_edna[,3:4],grps)

# Plot sites by kmeans clusters
cluster_coords$Cluster <- as.factor(cluster_coords$Cluster)
clustermap <- map2 + geom_point(data=cluster_coords,
                                aes(x=long,y=lat,
                                    col=Cluster),
                                size=3, inherit.aes=F)+
  scale_color_manual(values=c("darkblue","orange"))+
  theme(legend.position="right") +
  theme(axis.title=element_text(size=18))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=18))+ #change leg
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))
#ggtitle("Total number of taxa (species-level) detected per site per year")
clustermap



#---------- Additional figures---------
# Map of study sites coloured by year
cert3_geog <- read.csv(here("R output files", "cert3_geog.csv"))
sites_tomap <- cert3_geog %>% group_by(site,year) %>% distinct(Latitude,Longitude)
sites_tomap$year <- as.factor(sites_tomap$year)
sitemap <- map2 + geom_point(data=sites_tomap,
                             aes(x=Longitude,y=Latitude,
                                 col=year),
                             size=3, alpha=0.7,inherit.aes=F)+
  scale_colour_manual(values=c("darkred","darkgreen","darkblue"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.title = element_text(size = 14,colour="black"))+
  theme(axis.text=element_text(size=12, colour="black"))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14))+ #change legend text font size
  theme(legend.position="right")
sitemap

# Map of study sites (coloured by fish region)
fishregs_tomap <- edna4 %>% group_by(lat,long) %>% distinct(fish_reg)
colnames(fishregs_tomap)[3] <- "Region"

fishregmap <- map2 + geom_point(data=fishregs_tomap,
                                aes(x=long,y=lat,
                                    col=Region),
                                size=3, alpha=0.7,inherit.aes=F)+
  scale_color_manual(values=c("Newfoundland"="maroon",
                              "Southeast Labrador"="darkorange",
                              "Lake Melville"="red",
                              "Northern Labrador (to Fraser River)" ="blue",
                              "Northern Labrador (above Fraser River)"="darkgreen"))+
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme(axis.title = element_text(size = 14,colour="black"))+
  theme(axis.text=element_text(size=12, colour="black"))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14))+ #change legend text font size
  theme(legend.position="right")
fishregmap


# Inset map of Canada:
library(maps)
library(mapdata)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)


ne_countries(scale=10,returnclass="sf") %>%
  filter(geounit=="Canada") %>%
  ggplot() +
  geom_sf(fill="white") +
  theme(panel.border = element_rect(color = "black",
                                    fill = NA,
                                    size = 1))+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )




