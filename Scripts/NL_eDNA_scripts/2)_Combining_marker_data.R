###########################################################################
#------Script for combining presence/absence data from all 3 markers -----------
# Also contains code for generating summary stats on the data (i.e. means depths, etc.)

###########################################################################


#--------- Combining dfs with depth/taxon/sample AND total depth/sample----
library(here)
teleo <- read.csv(here("R output files", "12Steleo_depths_perfieldsample.csv"))
mifish <- read.csv(here("R output files", "mifish_depths_perfieldsample.csv"))
fishe <- read.csv(here("R output files","fishe_depths_perfieldsample.csv"))

# Bind all 3 dfs together:
depth_all <- rbind(teleo,mifish,fishe)

# Convert depths for each taxon to a proportion of total depths
# Note: this is the corrected depth for each taxon (after subtracting reads from negs)
# out of total depth per marker (uncorrected, for ALL taxa (not just fish) that were hit with BLAST)
# Note that for 12S markers this was everything; but FISHE CO1 was limited to vertebrates
# 12S markers did not hit anything that was not a vertebrate so this should be comparable
depth_all$tax_prop <- depth_all$depth/depth_all$tot_depth

#Testing some thresholds:
#Drake et al. 2021 paper states optimal thresholds combined with "max contam" methods were
# 0.002- 0.005
#To test, first select only hits:
library(dplyr)
depth_all %>% filter(tax_prop>0) # Result: 5331 out of 50528 rows were detections
#Next: how many of these detections were >=0.002 of the total sample depth?
depth_all %>% filter(tax_prop>=0.002)
# 4820/5331 were above this threshold - seems pretty reasonable as a threshold
# This was lowest threshold that performed well in Drake et al. paper (appropriate since previous cleaning steps were quite conservative)
# Also take a look at what is not meeting this threshold:
depth_all %>% filter(tax_prop>0) %>% filter(tax_prop<0.002) %>% 
  filter(marker=="12Steleo")
# Getting rid of 338 12Steleo hits, 160 Mifish hits, 13 FISHE hits

# ----- Calculate proportion of samples/marker/taxon that are above the threshold----
# Distinguish between all hits (regardless of depth), and higher certainty (above the threshold)
# 1. Above the threshold:
# Load in year data and append year to df
library(dplyr)
year <- read.csv(here("year_site.csv"))
agg_findyear_full <- year %>% 
  group_by(site,year) %>% summarise_each(funs(sum)) %>% as.data.frame()

# Create stub names for unique sites
depth_all$stub <- substr(depth_all$site,1,nchar(depth_all$site)-1)
# this next part is to deal with the fact that some sample names have additional numbers at the end
# Need to remove an additional letter, but only from sample names longer than 3 letters
fourplus_depth <- depth_all %>% filter(nchar(stub)>3) #select sample names > 3 letters
three_depth <- depth_all %>% filter(nchar(stub) <=3) #select sample names= 3 letters
fourplus_depth$stub <- sub("[A-Z]$","",fourplus_depth$stub) #remove last character from sample names > 3 letters ONLY if a letter (Not number)
depth_all1 <- rbind(fourplus_depth,three_depth) #bind back together
colnames(depth_all1)[1] <- "sample"
colnames(depth_all1)[7] <- "site"

depth_year <- merge(depth_all1,agg_findyear_full,
                    by="site",all.x=TRUE,all.y=FALSE)
#Check which are missing year and add in manually:
depth_year %>% filter(is.na(year)) #OLL was coded as OL in cega data
depth_oll <- depth_year %>% filter(site=="OLL") %>% mutate(year=c("2020"))
depth_notoll <- depth_year %>% filter(site!="OLL")
depth_year <- rbind(depth_notoll,depth_oll)

#Convert NAs in depth_year$tax_prop to 0 (from dividing 0 depth/0 total depth)
depth_year[,7][is.na(depth_year[,7])] <- 0

#Create column for number of detections (per marker per site) that do/do not meet threshold:
cert <- depth_year %>%
  group_by(site,taxon,year,marker) %>%
  mutate(n_meet_thresh=sum(tax_prop>=0.002),
         n_not_thresh=sum(tax_prop>0 & tax_prop<0.002)) %>% 
  ungroup() %>% as.data.frame() 

# Check that this worked (shouldn't be any rows that have tax_prop >=2 and n_meet_thresh as NA)
cert %>% filter(tax_prop>=0.002) %>% filter(is.na(n_meet_thresh)) # All good yay

# Create column for total number of marker x sample combos per site:
# Grouping by site, taxon, year (but not marker) will give number of samples/site
# REGARDLESS of whether a marker was picked up in all samples or taxon was picked up in all samples
cert$sample <- as.factor(cert$sample)
cert1 <- cert %>%
  group_by(site,year) %>%
  mutate(n_tot=n_distinct(sample)) %>% ungroup() %>% as.data.frame()

# Add columns for proportion of samples/site/marker that meet cutoff
# and also for detections that don't:
cert1$prop_meet_thresh <- cert1$n_meet_thresh/cert1$n_tot
cert1$prop_not_thresh <- cert1$n_not_thresh/cert1$n_tot
write.csv(cert1,"cert1.csv",row.names=FALSE) # save for using in making circlize plots

# Find a consensus for each taxon/marker/site
cert2 <- cert1 %>% group_by(site,taxon,year,marker) %>% 
  summarise(prop_marker_consensus= list(unique(prop_meet_thresh)),
            prop_notthresh_consensus= list(unique(prop_not_thresh)),
            marker_taxon_depth=sum(depth)) %>%
  as.data.frame()

# Check for any prop_marker_consensus that have two values assigned:
# prop=0: n=14107
# prop=1: n=928
# prop=1/3: n=503
# prop=2/3: n=498
# prop=1/6: n=33
# prop=3/6: n=88
# prop=5/6: n=5
# prop=1/4: n=29
# prop=3/4: n=32
# prop=1/5: n=0
# prop=2/5: n=0
# prop=3/5: n=0
# prop=4/5: n=0
# total: n=16223

# Average proportion of samples/site that meet threshold across the 3 markers:
cert2$prop_marker_consensus <- unlist(cert2$prop_marker_consensus)
cert2$prop_markernotthresh_consensus <- unlist(cert2$prop_notthresh_consensus)

cert3 <- cert2 %>% group_by(site,taxon,year) %>%
  summarise(prop_site_consensus= mean(prop_marker_consensus),
            prop_site_notthresh_consensus=mean(prop_markernotthresh_consensus)) %>%
  as.data.frame()
# This averaging takes into account that some markers do not detect certain species
# Only averages across markers that are able to detect that taxon (i.e. is somewhere in dataset)

# Note: Myoxocephalus and Ammodytes picked up by FISHE and MIFISHU, 
# but nothing at species level - should include these in addition to species-level
# assignments in downstream analyses
# Everthing else resolved to genus or family has something that could correspond in 
# more specific taxonomy

#------Check for missing data and fix-------
# Create column for river code:
cert3$RiverCode <- substr(cert3$site, 0, 3)
# merge with site coordinate data:
edna_coords <- read.csv(here("../","Working Data", "site_coords_allyears.csv"))
colnames(edna_coords)[4] <- "year"
cert3_geog <- merge(cert3,edna_coords,by=c('RiverCode',
                                           'year'),all.x=TRUE,all.y=TRUE)
# Check for missing data:
cert3_geog %>% filter(is.na(Latitude))
cert3_geog %>% filter(is.na(taxon))
# CRK entered as CKR in coordinate data
# Fraser River has no coordinate data (did not have centre coords so was missed)
# First add in RiverCode1 column to sort by
cert3_geog$RiverCode1 <- substr(cert3_geog$RiverCode, 0, 3)
cert3_geog[cert3_geog$RiverCode1=="CRK", "Latitude"] <- 53.804400
cert3_geog[cert3_geog$RiverCode1=="CRK", "Longitude"] <- -60.839270
cert3_geog[cert3_geog$RiverCode1=="FRA", "Latitude"] <- 56.625720
cert3_geog[cert3_geog$RiverCode1=="FRA", "Longitude"] <- -62.371610
#Get rid of RiverCode1 column again
cert3_geog<- cert3_geog[,-10]
# There are 7 rows with coordinate data but no taxa- can ignore these since they won't be plotted anyway
# These 7 rows are comprised of the remaining CRK row, as well as 3 sites that had an
# additional 2 subsites (which were grouped all together under main site record)

write.csv(cert3_geog,"cert3_geog.csv",row.names=FALSE)

#-------------Creating supplementary table of all taxa occurrences by site and year----------
# Supplemental Table 3 in manuscript
library(here)
library(dplyr)
library(tidyr)
cert3_geog <- read.csv(here("R output files","cert3_geog.csv"))
# First, exclude sites with no taxon data (were grouped with their main site)
cert3_geog <- cert3_geog %>% filter(!is.na(taxon))
# select only high cert detections, convert to pres/abs (0/1) and create wide df:
dets_wide <- cert3_geog %>% select(RiverCode,year,site,taxon,prop_site_consensus,
                                   RiverName, Latitude,Longitude) %>%
  mutate(prop_site_consensus= case_when(prop_site_consensus==0 ~ " ",
                                        prop_site_consensus>0 ~ "X")) %>%
  group_by(RiverCode,year,site,taxon,prop_site_consensus,
           RiverName, Latitude,Longitude) %>% 
  spread(taxon, prop_site_consensus) %>% as.data.frame()
# Convert NAs to zero
dets_wide[,c(7:47)][is.na(dets_wide)[,c(7:47)]] <- " "
dets_wide <- dets_wide[,-c(3)]
dets_wide <- dets_wide %>% relocate(RiverName, .before= RiverCode)
write.csv(dets_wide, "eDNA_taxon_occurrence.csv",row.names=FALSE)

# -------Check for consistency in detections in Labrador sites sampled in both years-------
year2019 <- cert3_geog %>% filter(year=="2019") %>% distinct(RiverCode)
year2021 <- cert3_geog %>% filter(year=="2021") %>% distinct(RiverCode)
both <- year2019 %>% filter(RiverCode %in% year2021$RiverCode)

# Which species were detected at the sites common to 2019 and 2021?
dets2019 <- cert3_geog %>% filter(RiverCode=="109" | RiverCode=="CHR" | RiverCode=="ENG" |
                                    RiverCode=="FOR" | RiverCode=="HUR" | RiverCode=="IKA" |
                                    RiverCode=="KIN" | RiverCode=="MBB" | RiverCode=="PAN" |
                                    RiverCode=="PIN" | RiverCode=="SAN" | RiverCode=="SUS") %>%
  filter(prop_site_consensus>0) %>% filter(year=="2019") %>% distinct(taxon)

dets2021 <- cert3_geog %>% filter(RiverCode=="109" | RiverCode=="CHR" | RiverCode=="ENG" |
                                    RiverCode=="FOR" | RiverCode=="HUR" | RiverCode=="IKA" |
                                    RiverCode=="KIN" | RiverCode=="MBB" | RiverCode=="PAN" |
                                    RiverCode=="PIN" | RiverCode=="SAN" | RiverCode=="SUS") %>%
  filter(prop_site_consensus>0) %>% filter(year=="2021") %>% distinct(taxon)

bothtaxa <- dets2019 %>% filter(taxon %in% dets2021$taxon)
notboth <- dets2019 %>% filter(! taxon %in% dets2021$taxon)

#---------------General data summary-------------
# read depth per sapmle (both clean (i.e. post-bioinformatics), and raw)
fishdepth <- depth_year %>% group_by(sample,marker) %>% 
  summarise(clean_depth=sum(depth))
fishdepth %>% group_by(marker) %>% summarise(mean_depth=mean(clean_depth))
# mean depths (fish asvs, after bioinformatic filtering including blank subtraction)
# 12Steleo: 669940
# MIFISHU: 937463
# FISHE: 4857 

rawfishdepth <- depth_year %>% group_by(sample,marker) %>% 
  summarise(tot_depth=mean(tot_depth))
rawfishdepth %>% group_by(marker) %>% summarise(mean_depth=mean(tot_depth))
# mean raw depths (all asvs)
# 12Steleo: 881837
# MIFISHU: 964985
# FISHE: 5779 (note this is only vertebrate ASVs still!)

# Numbers of ASVs per sample and marker
# Start out with raw ASVs
# 1. 12Steleo marker
library(here)
# a) number of raw asvs (before blasting)
asv_12s_full <- read.csv(here("Compute Canada Transfers",
                              "12Steleo","asv_counts_12steleo_full.csv"))
# 12Steleo marker had 7868 ASVs

# b) number of ASVs that assigned at at least 99% sss
teleo_asv99 <- read.csv(here("R output files",
                             "field_samples_12steleo.csv"))
# 57 total asvs assigned at at least 99%

# c) Only fish ASVs, cleaned 
teleo_fish_asvs <- read.csv(here("R output files",
                                 "teleo_asvs_persample.csv"))
onlyfish_teleo_asvs <- teleo_fish_asvs[,c(1,2,3,4,6,8,10,12,13,14,18,21,23,24,
                                          26:29,33,35,36,38,
                                          41:43,47,51:53,57)]
# total of 29 fish ASVs for 12Steleo marker
onlyfish_teleo_asvs[is.na(onlyfish_teleo_asvs)] <- 0 #convert NAs to zeros
onlyfish_teleo_asvs <- data.frame(rowSums(onlyfish_teleo_asvs[,2:30]>0))
colnames(onlyfish_teleo_asvs)[1] <- "teleo_fish_asvs"   
mean(onlyfish_teleo_asvs$teleo_fish_asvs)
min(onlyfish_teleo_asvs$teleo_fish_asvs)
max(onlyfish_teleo_asvs$teleo_fish_asvs)
# mean number ASVs/sample for FISHE marker is 3.881 (range 0-10)

# 2. MIFISHU marker
# a) number of raw asvs (before blasting)
asv_mifish_full <- read.csv(here("Compute Canada Transfers",
                                 "MIFISHU","asv_counts_mifishu_full.csv"))
# MIFISHU marker had 670 ASVs

# b) number of ASVs that assigned at at least 99%
mifish_asv99 <- read.csv(here("R output files",
                              "field_samples_mifishu.csv"))
# 101 total asvs assigned at at least 99%

# c) Only fish ASVs, cleaned
mifish_fish_asvs <- read.csv(here("R output files",
                                  "mifish_asvs_persample.csv"))
onlyfish_mifish_asvs <- mifish_fish_asvs[,c(1,2,3,6,8:10,15,16,21,23,
                                            26,28,30,35,36,40:42,44,46,
                                            48:50,54,55,58,59,62,63,
                                            70,74,79,80,90,91,97)]
# total of 36 fish ASVs for MIFISHU marker
onlyfish_mifish_asvs[is.na(onlyfish_mifish_asvs)] <- 0 #convert NAs to zeros
onlyfish_mifish_asvs <- data.frame(rowSums(onlyfish_mifish_asvs[,2:37]>0))
colnames(onlyfish_mifish_asvs)[1] <- "mifish_fish_asvs" 
mean(onlyfish_mifish_asvs$mifish_fish_asvs)
min(onlyfish_mifish_asvs$mifish_fish_asvs)
max(onlyfish_mifish_asvs$mifish_fish_asvs)
# mean number ASVs/sample for FISHE marker is 4.26 (range 0-11)

# 3. FISHE marker
# a) number of raw asvs (before blasting)
asv_fishe_full <- read.csv(here("Compute Canada Transfers",
                                "FISHE","reduced_asv_counts.csv"))
# FISHE marker had 751 total ASVs

# b) number of asvs that assigned at at least 99%
fishe_asv99 <- read.csv(here("R output files",
                             "field_samples_fishe.csv"))
# 121 total asvs assiged at at least 99%

# c) only fish asvs, cleaned
fishe_fish_asvs <- read.csv(here("R output files",
                                 "fishe_asvs_persample.csv"))
onlyfish_fish_asvs <- fishe_fish_asvs[,c(1,12,13,17,21,24,26,28:31,
                                         33,35,36,38,39,42,44,45,50:52,
                                         55:57,61,73,75,78,80,82,84,
                                         87:90,92,94,98,100,103,
                                         106,110,111,121)]

# total of 44 fish ASVs for FISHE marker
onlyfish_fish_asvs[is.na(onlyfish_fish_asvs)] <- 0 #convert NAs to zeros
onlyfish_fish_asvs <- data.frame(rowSums(onlyfish_fish_asvs[,2:45]>0))
colnames(onlyfish_fish_asvs)[1] <- "fishe_fish_asvs"   
mean(onlyfish_fish_asvs$fishe_fish_asvs)
min(onlyfish_fish_asvs$fishe_fish_asvs)
max(onlyfish_fish_asvs$fishe_fish_asvs)
# mean number ASVs/sample for FISHE marker is 2.875 (range 0-12)

# Mean depth per species per marker
depth_species_marker <- cert1 %>% group_by(marker,taxon) %>% filter(depth>0) %>%
  filter(grepl(" ",taxon)) %>%
  summarise(mean_depth=mean(depth),stdev_depth=sd(depth))
depth_species_marker$mean_depth <- round(depth_species_marker$mean_depth,0)
depth_species_marker$stdev_depth <- round(depth_species_marker$stdev_depth,0)
depth_species_marker <- depth_species_marker %>% group_by(marker) %>%
  arrange(desc(mean_depth),.by_group = TRUE)

# Figure 2 panel in manuscript:
ggplot(depth_species_marker, aes(x = reorder(taxon, mean_depth),
                                 y = mean_depth, color = marker)) +
  geom_point(size=4) +
  scale_colour_manual(values = c("darkblue", "cornflowerblue", "cyan4"))+
  geom_errorbar(aes(ymax = (mean_depth + stdev_depth), 
                    ymin = mean_depth, color=marker), size = .5, width=0.3) +
  scale_y_continuous(trans="log10") +
  #scale_y_continuous(trans=scales::pseudo_log_trans(base = 10))+
  theme_bw() + 
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank() )+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size=18,colour="black"))+
  theme(axis.title=element_text(size=18))+
  theme(axis.text.y = element_text(size=18,colour="black"))+
  theme(axis.line.x = element_line(color="black", size = 1),
        axis.line.y = element_line(color="black", size = 1))+
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=18), #change legend title font size
        legend.text = element_text(size=18))+ #change legend text font size
  labs(y="Mean depth", x="Taxon")




