###########################################################################
# ------Script for environmental association analyses of eDNA community composition data,
# -------compared with historical data (PADA)----

###########################################################################


# -------------PADA vs. eDNA comparisons------------
# load in eDNA data (edna3.csv, from edna_comm_comp.R)
library(here)
library(dplyr)
edna3 <- read.csv(here("R output files","edna3.csv"))
# load in pada data (river sites only)
pada_river2 <- read.csv(here("R output files","pada_river.csv"))

# Match each entry in pada_river2 to nearest site in edna3
library(hutilscpp)
river_match <- match_nrst_haversine(pada_river2$lat, 
                                    pada_river2$long,
                                    edna3$lat,
                                    edna3$long,
                                    close_enough="0.25km",
                                    .verify_box = TRUE)

# bind back with original pada_river df (in same order as started):
river_match <- cbind(pada_river2, river_match)
# Bind with matching entry in edna3 df
edna3$pos <- as.integer(row.names(edna3))
river_merge <- merge(river_match, edna3, by="pos")

# Select only those with distance < 25km
river_merge25 <- river_merge %>% filter(dist<=25)

# Keep only needed columns
river_merge25b <- river_merge25 %>%
  select(location.x, year.x,lat.x,long.x,spec_name.x,ecoreg.y,subreg.y,
         source.x,fish_reg.y,location.y,lat.y,long.y) %>%
  rename(pada_location=location.x,pada_year=year.x,pada_lat=lat.x,
         pada_long=long.x,pada_spec_name=spec_name.x,
         ecoreg=ecoreg.y,subreg=subreg.y,source=source.x,
         fish_reg=fish_reg.y,edna_location=location.y,
         edna_lat=lat.y,edna_long=long.y) 
# note that ecoreg, subreg, and fish_reg chosen are for eDNA coords

# Data summary:
river_merge25b %>% filter(pada_year < 1987) %>%
  distinct(edna_location) 
# This keeps 104 "paired" locations that are river sites only for pre-1987
obs <- river_merge25b %>% filter(pada_year <1987) %>%
  group_by(edna_lat, edna_long,edna_location) %>%
  tally() 
obs %>% filter(n<5) # 44 locations have <5 observations 
obs %>% filter(n>=5 & n<10) # 15 locations have between 5 and 10 observations
obs %>% filter(n>=10 & n<50) # 22 locations have between 10 and 50 observations
obs %>% filter(n>50) # 27 locations have >50 observations (some have thousands) 

# Make a new df with just pre-1987 data to compare with eDNA:
pre87_river <- river_merge25b %>%
  filter(pada_year<1987) %>%
  select(edna_lat,edna_long,
         pada_spec_name,source,fish_reg,edna_location) %>%
  rename(lat=edna_lat,
         long=edna_long,spec_name=pada_spec_name,RiverCode=edna_location) %>%
  mutate(source="historical")

# Choose same edna columns and bind together with pre87_river
edna4 <- edna3[,-10] #get rid of pos column now
edna4 <- edna4 %>% select(location,lat,long,
                          spec_name,source,fish_reg) %>%
  rename(RiverCode=location)
river_comp <- rbind(pre87_river,edna4)

# Convert to wide format:
library(tidyr)
river_comp_wide <- river_comp %>%
  select(lat, long, RiverCode,spec_name,source,fish_reg) %>%
  group_by(lat, long,RiverCode,spec_name,source,fish_reg) %>%
  summarise(n_obs=n()) %>%
  spread(spec_name,n_obs) %>%
  as.data.frame() %>%
  replace(is.na(.),0)
# Convert observation counts to just 1 (pres/abs format)
river_comp_wide <- as.data.frame(cbind(river_comp_wide[,1:5],(+(river_comp_wide[-(1:5)] > 0))))

# Split dfs apart again:
edna_wide <- river_comp_wide %>% filter(source=="eDNA")
pada87_wide <- river_comp_wide %>% filter(source=="historical")

# Keep only unique lat/longs in edna_wide that also exist in pada87_wide
edna_final <- subset(edna_wide,
                     lat %in% unique(pada87_wide$lat))
rownames(edna_final) <- rep(1:108)

#get rid of last column in each ("Unspecified stickleback")
edna_final <- edna_final %>% 
  select(-c(`unspecified stickleback`))
pada87_final <- pada87_wide %>% 
  select(-c(`unspecified stickleback`))
# Check for any sites that may have lost their only detection through exclusions above:
rowSums(edna_final[,6:38]) # no zero-sum rows
rowSums(pada87_final[,6:38]) # no zero-sum rows
rownames(edna_final) <- rep(1:108)
rownames(pada87_final) <- rep(1:108)


# Arrange both dfs by lat so they are in the same order:
edna_final <- edna_final %>% arrange(lat)
pada87_final <- pada87_final %>% arrange(lat)
write.csv(pada87_final,"pada87_final.csv",row.names=FALSE) # save for using to generate MEMs in MEMs.R
write.csv(edna_final,"edna_final_comp.csv", row.names=FALSE)

# -----------Temporal beta diversity analysis of all "paired" sites--------
library(adespatial)
all_tbi <- TBI(pada87_final[,c(6:38)],edna_final[,c(6:38)],
               method="jaccard",   #jaccard dissimilarity index
               pa.tr=TRUE,         #data transformed to pres/abs (already there but do it anyway?) 
               nperm=9999,         #number of permutations
               BCD=TRUE,           #compute and present B and C components
               replace=TRUE,       # sampling with replacement (bootstrapping)
               test.BC=TRUE,       #t-test of species gains and losses
               test.t.perm=TRUE,   #permutation test of gains and losses 
               seed=1234,          #seed set so reproducible
               clock=TRUE)  

# --------Plot of losses and gains -----
library(ggplot2)
change_df <- cbind(pada87_final[,c(3,5)],all_tbi$BCD.mat)
colnames(change_df) <- c("RiverCode","Region","Losses","Gains",
                         "Difference","Change")
ggplot(change_df,aes(x=Losses,y=Gains,color=Region))+
  geom_point(size=4)+
  geom_hline(yintercept=0,linetype="dotted")+
  geom_vline(xintercept=0,linetype="dotted") +
  coord_fixed()+
  geom_abline(intercept=0, slope=1)+
  scale_color_manual(values=c("Newfoundland"="maroon",
                              "Southeast Labrador"="darkorange",
                              "Lake Melville"="red",
                              "Northern Labrador (to Fraser River)" ="blue",
                              "Northern Labrador (above Fraser River)"="darkgreen"))+
  theme_bw()+
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(axis.title = element_text(size = 14,colour="black"))+
  theme(axis.text=element_text(size=12, colour="black"))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=16), #change legend title font size
        legend.text = element_text(size=14)) #change legend text font size
#facet_wrap(~fish_reg)+
#ggtitle("Species Losses and Gains (historical vs. eDNA)")

#---------Generating MEMs with marine distances for historical subset of data------
pada87_final <- read.csv(here("R output files","pada87_final.csv"))
coords_pada <- pada87_final %>% select(lat,long)

bathydata2 <- getNOAA.bathy(-68,-50,46,62, res=2,keep=T)
trans1 <- marmap::trans.mat(bathydata2)


# Compute distance from each site to nearest point at the coast
# isobath=0 is nearest point on coast, so use isobath= -10
mar_dist_pada <- dist2isobath(bathydata2, coords_pada$long, coords_pada$lat, isobath = -10)
# Check depths of new coordinates:
onland_pada <- get.depth(bathydata2, mar_dist_pada[,4:5], locator = FALSE) %>%
  filter(depth>=0) # 13 points still have a positive depth:

notland_pada <- get.depth(bathydata2, mar_dist_pada[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start_pada <- mar_dist_pada %>% filter(end.lat %in% notland_pada$lat)

# run dist2isobath again for these 28 points:
land_pada <- mar_dist_pada %>% filter(end.lat %in% onland_pada$lat) # filter for starting coords that were still not on land
mar_dist_pada2 <- dist2isobath(bathydata2, land_pada$start.lon, land_pada$start.lat, isobath = -20)
onland_pada2 <- get.depth(bathydata2, mar_dist_pada2[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # 5 of them are still on land
# add sites that are not not on land to notland df:
notland_pada2 <- get.depth(bathydata2, mar_dist_pada2[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start_pada2 <- mar_dist_pada2 %>% filter(end.lat %in% notland_pada2$lat)
notland_start_pada <- rbind(notland_start_pada,notland_start_pada2)

# Bump up isobath again:
land_pada2 <- mar_dist_pada2 %>% filter(end.lat %in% onland_pada2$lat) # filter for starting coords that were still not on land
mar_dist_pada3 <- dist2isobath(bathydata2, land_pada2$start.lon, land_pada2$start.lat, isobath = -30)
onland_pada3 <- get.depth(bathydata2, mar_dist_pada3[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # 3 still on land, bump up
# add sites that are not not on land to notland df:
notland_pada3 <- get.depth(bathydata2, mar_dist_pada3[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start_pada3 <- mar_dist_pada3 %>% filter(end.lat %in% notland_pada3$lat)
notland_start_pada <- rbind(notland_start_pada,notland_start_pada3)

# Bump up isobath again:
land_pada3 <- mar_dist_pada3 %>% filter(end.lat %in% onland_pada3$lat) # filter for starting coords that were still not on land
mar_dist_pada4 <- dist2isobath(bathydata2, land_pada3$start.lon, land_pada3$start.lat, isobath = -50)
onland_pada4 <- get.depth(bathydata2, mar_dist_pada4[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # 2 still on land
# add sites that are not not on land to notland df:
notland_pada4 <- get.depth(bathydata2, mar_dist_pada4[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start_pada4 <- mar_dist_pada4 %>% filter(end.lat %in% notland_pada4$lat)
notland_start_pada <- rbind(notland_start_pada,notland_start_pada4)

# Bump up isobath again:
land_pada4 <- mar_dist_pada4 %>% filter(end.lat %in% onland_pada4$lat) # filter for starting coords that were still not on land
mar_dist_pada5 <- dist2isobath(bathydata2, land_pada4$start.lon, land_pada4$start.lat, isobath = -70)
land_pada5 <- get.depth(bathydata2, mar_dist_pada5[,4:5], locator = FALSE) %>%
  dplyr::filter(depth>=0) # None are on land anymore
# add sites that are not not on land to notland df:
notland_pada5 <- get.depth(bathydata2, mar_dist_pada5[,4:5], locator = FALSE) %>%
  filter(depth<0)
notland_start_pada5 <- mar_dist_pada5 %>% filter(end.lat %in% notland_pada5$lat)
notland_start_pada <- rbind(notland_start_pada,notland_start_pada5)
# notland_start_pada contains the starting site coordinates, as well as 
# the coordinates that are the closest in-ocean coordinates

notland_start_pada <- notland_start_pada %>% arrange(start.lat) #arrange by starting latitude to match up with community data later
marine_sites_pada <- notland_start_pada[,4:5] 
lcp_pada <- lc.dist(trans1, marine_sites_pada, res = "dist", meters=T)

coords_pada <- coords_pada %>% arrange(lat)
nb_pada = spdep::graph2nb(spdep::gabrielneigh(coords_pada), sym = TRUE)  # Gabriel graph: neighbor definition (a list of neighbors for each site/individual)
spdep::plot.nb(nb_pada, coords_pada, pch = 21)
disttri_pada = spdep::nbdists(nb_pada, coords_pada,longlat=TRUE)  # Add longlat=T for lat/long coordinates
# modify the weights based on the actual geographic (marine) distances between neighbors
im_pada = 1/lcp_pada  # Use inverse distance weights, (the further away, the lower the value)
im_pada[!is.finite(im_pada)]<- 0
im_pada2 <- as.matrix(im_pada)

# get a revised spatial weights matrix
fdist_pada = lapply(1:length(nb_pada), function(x){c(sapply(nb_pada[x], 
                                                            function(z){ im_pada2[x,z]}))}) # end lapply
listW_pada = spdep::nb2listw(nb_pada, style="W", glist = fdist_pada) #glist adds revised weights
summary(sapply(listW_pada$weights, sum))
mems_pada = adespatial::scores.listw(listW_pada, MEM.autocor = "positive") 
summary(mems_pada)
# Generates 49 positive MEMs

# merge community data with MEMs
# MEMs were arranged by latitude, so start by arranging each community df by latitude too
pada87_final <- pada87_final %>% arrange(lat)
pada87_mem <- cbind(pada87_final,mems_pada)
edna_final <- read.csv(here("R output files","edna_final_comp.csv"))
edna_final <- edna_final %>% arrange(lat)
edna_mem <- cbind(edna_final,mems_pada)
full_comm_marinemem <- rbind(pada87_mem,edna_mem) # merge eDNA and pada dfs together
#merge with climate data:
all_climate2 <- read.csv(here("R output files","climate_avgs.csv"))
comm_env_comp_marinemem <- merge(full_comm_marinemem, all_climate2,by=c("RiverCode","source"),
                                 all.x=TRUE,all.y=FALSE)
# Save for use in actual RDAs too
write.csv(comm_env_comp_marinemem,"comm_env_comp_marinemem.csv",row.names=FALSE)


# Check for correlations among explanatory variables (MEMs and climate)
cors_marinemem <- data.frame(cor(comm_env_comp_marinemem[,39:94]))
#First check if any MEMs are strongly correlated (>0.75) with climate vars:
cors_clim_marinemem <- cors_marinemem[1:49,50:56]
cors_clim_marinemem %>% filter_all(any_vars(abs(.) > abs(0.75)))
# There are no MEMs correlated with climate variables above 0.75

# Run RDA wit just MEMs and use ordistep to pick a subset for use in RDA
full_rda_mems_pada <- comm_env_comp_marinemem[,39:87]
library(vegan)
dbrda_edna_mems_pada <- capscale(comm_env_comp_marinemem[,c(6:38)] ~ .,
                                 full_rda_mems_pada,
                                 dist="jaccard")
summary(dbrda_edna_mems_pada)
# this rda explains 44.37% of the variation in the data
#https://www.nature.com/articles/s41437-020-0352-6

# Compare with null rda, which is just a pcoa
dbrda_edna_null_pada <- capscale(comm_env_comp_marinemem[,c(6:38)] ~ 1,
                                 full_rda_mems_pada,
                                 dist="jaccard")
selection_pada <- ordiR2step(dbrda_edna_null_pada, scope = formula(dbrda_edna_mems_pada), 
                             direction="both")
selection_pada$anova  # this gives significant MEMs, 17 in total
# save selected MEMs 
marine_edna_mems_pada <- data.frame(cbind(full_rda_mems_pada$MEM6,full_rda_mems_pada$MEM5,
                                          full_rda_mems_pada$MEM4,
                                          full_rda_mems_pada$MEM8,full_rda_mems_pada$MEM3,
                                          full_rda_mems_pada$MEM7,
                                          full_rda_mems_pada$MEM12,full_rda_mems_pada$MEM46,
                                          full_rda_mems_pada$MEM9,
                                          full_rda_mems_pada$MEM31,full_rda_mems_pada$MEM28,
                                          full_rda_mems_pada$MEM38,
                                          full_rda_mems_pada$MEM1,full_rda_mems_pada$MEM14,
                                          full_rda_mems_pada$MEM16,
                                          full_rda_mems_pada$MEM40,full_rda_mems_pada$MEM19))
colnames(marine_edna_mems_pada) <- c("MEM6","MEM5",
                                     "MEM4",
                                     "MEM8","MEM3",
                                     "MEM7",
                                     "MEM12","MEM46",
                                     "MEM9",
                                     "MEM31","MEM28",
                                     "MEM38",
                                     "MEM1","MEM14",
                                     "MEM16",
                                     "MEM40","MEM19")
write.csv(marine_edna_mems_pada,"marine_edna_mems_pada.csv",row.names=FALSE) # MEMs are still arranged by latitude



#-------------dbRDAs using climate variables, data source, and MEMs corrected for marine distances------------------
# Load in sorted community and climate df for comparison, and selected MEMs
library(here)
library(dplyr)
comm_env_comp_marinemem <- read.csv(here("R output files","comm_env_comp_marinemem.csv"))
marine_edna_mems_pada <- read.csv(here("R output files","marine_edna_mems_pada.csv"))

# merge community data with MEMs
comm_env_comp_mar <- cbind(comm_env_comp_marinemem[,c(1:38,88:94)],marine_edna_mems_pada)
comm_env_comp_mar$source <- as.factor(comm_env_comp_mar$source)

# For variance partitioning, need to run multiple models:
# 1. Full model:
full_rda_expl_pada_mar <- comm_env_comp_mar[,c(2,39:62)]
dbrda_edna_full_pada_mar <- capscale(comm_env_comp_mar[,c(6:38)] ~ mean_max_air_temp +
                                       mean_rel_humid + mean_wind_speed_2m +
                                       mean_solar_rad + mean_atm_press +
                                       mean_snow_precip + mean_snow_accum +
                                       MEM6 + MEM5 +
                                       MEM4 +
                                       MEM8 + MEM3 +
                                       MEM7 +
                                       MEM12 + MEM46 +
                                       MEM9 +
                                       MEM31 + MEM28 +
                                       MEM38 +
                                       MEM1 + MEM14 +
                                       MEM16 +
                                       MEM40 + MEM19 + source, 
                                     full_rda_expl_pada_mar,
                                     dist="jaccard")
summary(dbrda_edna_full_pada_mar)
# total proportion of variance explained is 46.38% (unexplained is 53.62%)
# Inertia: total is 64.37, constrained is 29.85, unconstrained is 34.52
anova.cca(dbrda_edna_full_pada_mar)
# significant at p= 0.001, F=6.5726

# 2. Pure climate influence:
dbrda_edna_full_pada_pureclim <- capscale(comm_env_comp_mar[,c(6:38)] ~ mean_max_air_temp +
                                            mean_rel_humid + mean_wind_speed_2m +
                                            mean_solar_rad + mean_atm_press +
                                            mean_snow_precip + mean_snow_accum +
                                            Condition(source + MEM6 + MEM5 +
                                                        MEM4 +
                                                        MEM8 + MEM3 +
                                                        MEM7 +
                                                        MEM12 +MEM46 +
                                                        MEM9 +
                                                        MEM31+ MEM28 +
                                                        MEM38 +
                                                        MEM1 +MEM14 +
                                                        MEM16 +
                                                        MEM40 +MEM19), 
                                          full_rda_expl_pada_mar,
                                          dist="jaccard")

smry_dbrda_comp_pureclim <- summary(dbrda_edna_full_pada_pureclim)

# conditioned vars (17 significant MEMs + data source) explain 39.18%
# constrained vars (climate) explain 7.19%
anova.cca(dbrda_edna_full_pada_pureclim) # model is significant (p=0.001)
anova.cca(dbrda_edna_full_pada_pureclim,by="terms", permu=9999)
anova.cca(dbrda_edna_full_pada_pureclim,by="axis", permu=9999)
# all climate vars are significant at p<0.05
# first 3 axes are significant at p<0.05

# 2.5 Full Climate
# This one is not necessary for variance partitioning, but useful for visualization
# ie. it shows ALL climate influence (pure AND spatially-structured) AND source
dbrda_edna_full_pada_allclim <- capscale(comm_env_comp_mar[,c(6:38)] ~ source + mean_max_air_temp +
                                           mean_rel_humid + mean_wind_speed_2m +
                                           mean_solar_rad + mean_atm_press +
                                           mean_snow_precip + mean_snow_accum, 
                                         full_rda_expl_pada_mar,
                                         dist="jaccard")

smry_dbrda_comp_allclim <- summary(dbrda_edna_full_pada_allclim)

# 3. Pure space
dbrda_edna_full_pada_puregeo <- capscale(comm_env_comp_mar[,c(6:38)] ~  MEM6 + MEM5 +
                                           MEM4 +
                                           MEM8 + MEM3 +
                                           MEM7 +
                                           MEM12 + MEM46 +
                                           MEM9 +
                                           MEM31 + MEM28 +
                                           MEM38 +
                                           MEM1 + MEM14 +
                                           MEM16 +
                                           MEM40 + MEM19 +
                                           Condition(source + mean_max_air_temp +
                                                       mean_rel_humid + mean_wind_speed_2m +
                                                       mean_solar_rad + mean_atm_press +
                                                       mean_snow_precip + mean_snow_accum), 
                                         full_rda_expl_pada_mar,
                                         dist="jaccard") 

summary(dbrda_edna_full_pada_puregeo)
anova.cca(dbrda_edna_full_pada_puregeo)
# Conditioned variables (climate + data source) explain 36.61%
# Constraining variables (space) explain 9.77%

# 4. Pure data source:
dbrda_edna_full_pada_puresource <- capscale(comm_env_comp_mar[,c(6:38)] ~  source +
                                              Condition(mean_max_air_temp +
                                                          mean_rel_humid + mean_wind_speed_2m +
                                                          mean_solar_rad + mean_atm_press +
                                                          mean_snow_precip + mean_snow_accum +
                                                          MEM6 + MEM5 +
                                                          MEM4 +
                                                          MEM8 + MEM3 +
                                                          MEM7 +
                                                          MEM12 + MEM46 +
                                                          MEM9 +
                                                          MEM31 + MEM28 +
                                                          MEM38 +
                                                          MEM1 + MEM14 +
                                                          MEM16 +
                                                          MEM40 + MEM19), 
                                            full_rda_expl_pada_mar,
                                            dist="jaccard") 
summary(dbrda_edna_full_pada_puresource)
# conditioned variables (climate and MEMs) explain 41.98%
# constrained (data source) explains 4.40%
anova.cca(dbrda_edna_full_pada_puresource)

#---------Biplots for RDAs---------------
#save site loadings for first two RDA axes (second RDA):
site_loadings_rda_full_pada_mar2 <- as.data.frame(smry_dbrda_comp_pureclim$sites[,1:2])
vars_loadings_rda_full_pada_mar2 <- as.data.frame(smry_dbrda_comp_pureclim$biplot[,1:2])
rownames(vars_loadings_rda_full_pada_mar2) <- c("Max Air Temp",
                                                "Relative Humidity",
                                                "Wind Speed 2m",
                                                "Solar Radiation",
                                                "Atmospheric Pressure",
                                                "Snow Precip",
                                                "Snow Accum")
#bind back with comm_env_pc
climate_meta_rda_full_pada_mar2 <- cbind(comm_env_comp_mar[,1:5],
                                         site_loadings_rda_full_pada_mar2)

climate_meta_rda_full_pada_mar2 <- climate_meta_rda_full_pada_mar2 %>% 
  arrange(source, RiverCode)


# Biplot
library(ggplot2)
colnames(climate_meta_rda_full_pada_mar2)[5] <- "Region"
ggplot(climate_meta_rda_full_pada_mar2,aes(x=CAP1,y=CAP2, col=Region))+
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
  # coord_fixed() +
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
ggrepel::geom_label_repel(data=vars_loadings_rda_full_pada_mar2,aes(x=CAP1,y=CAP2,
                                                                    label=rownames(vars_loadings_rda_full_pada_mar2)),
                          # hjust=0.5*(1-sign(CAP1)),
                          # vjust=0.5*(1-sign(CAP2))),
                          color="black",size=3.5,alpha=0.8,position=position_dodge(0),
                          max.overlaps=Inf)+
  geom_segment(data=vars_loadings_rda_full_pada_mar2, 
               aes(x=0,xend=CAP1,y=0,yend=CAP2),
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
        legend.text = element_text(size=14))+ 
  xlab("CAP1") + ylab("CAP2")
#facet_wrap(~fish_reg)+
# ggtitle("db-RDA community composition (PADA vs. eDNA) and climate with MEMs corrected for marine distance")


#save site loadings for first two RDA axes (all enviro, RDA 2.5):
site_loadings_rda_nocond_mar <- as.data.frame(smry_dbrda_comp_allclim$sites[,1:2])
vars_loadings_rda_nocond_mar <- as.data.frame(smry_dbrda_comp_allclim$biplot[,1:2])
rownames(vars_loadings_rda_nocond_mar) <- c("Data Source","Max Air Temp",
                                            "Relative Humidity",
                                            "Wind Speed 2m",
                                            "Solar Radiation",
                                            "Atmospheric Pressure",
                                            "Snow Precip",
                                            "Snow Accum")
#bind back with comm_env_pc
climate_meta_rda_nocond_mar <- cbind(comm_env_comp_mar[,1:5],
                                     site_loadings_rda_nocond_mar)

climate_meta_rda_nocond_mar <- climate_meta_rda_nocond_mar %>% 
  arrange(source, RiverCode)


# Biplot
library(ggplot2)
colnames(climate_meta_rda_nocond_mar)[5] <- "Region"
ggplot(climate_meta_rda_nocond_mar,aes(x=CAP1,y=CAP2, col=Region))+
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
  # coord_fixed() +
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
ggrepel::geom_label_repel(data=vars_loadings_rda_nocond_mar,aes(x=CAP1,y=CAP2,
                                                                label=rownames(vars_loadings_rda_nocond_mar)),
                          #   hjust=0.5*(1-sign(CAP1)),
                          #   vjust=0.5*(1-sign(CAP2))),
                          color="black",size=3.5,alpha=0.8,position=position_dodge(0),
                          max.overlaps = Inf)+
  geom_segment(data=vars_loadings_rda_nocond_mar, 
               aes(x=0,xend=CAP1,y=0,yend=CAP2),
               color="black",size=1)+
  theme_bw()+
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
        legend.text = element_text(size=14))+ 
  xlab("CAP1") + ylab("CAP2")
#facet_wrap(~fish_reg)+
#ggtitle("db-RDA community composition (PADA vs. eDNA) and climate")


