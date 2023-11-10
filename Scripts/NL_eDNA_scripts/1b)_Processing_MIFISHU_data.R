
###########################################################################
#------Script for Processing BLAST output from MIFISHU marker data -----------

###########################################################################

# -----BLAST data processing-------------------
# Load in BLAST output data (transferred from Compute Canada)
library(here)
blastmf_full <- read.csv(here("Compute Canada Transfers", "MIFISHU", 
                              "blast_MIFISHU_full.csv"))

# Create column for "sequence similarity score (sss) based on percent identity (pident) and percent query coverage (qcovs)
blastmf_full$sss <- (blastmf_full$pident)*(blastmf_full$qcovs)/100

# Create new data frame only selecting rows with >=99% sss and number of different taxonomic hits ("counts)
library(dplyr)
speciesmf_full <- filter(blastmf_full,sss>=99)

# Create data frame for each ASV with sss>=99 and number of different taxonomic hits ("count")
species_per_read_full_mf <- data.frame(speciesmf_full %>% group_by(readID) %>%
                                         summarise(count = n_distinct(ssciname)))

# Create data frame listing each taxonomic hit for each ASV with sss>=99
speciesID_per_read_full_mf <- data.frame(speciesmf_full %>% group_by(readID) %>%
                                           distinct(ssciname))

#----- Assign BLAST ASV IDs to reads at each site-----
# Start easy; select readIDs with only one taxonomic match
one_match_mf <- data.frame(species_per_read_full_mf %>% filter(count==1))
# Merge this data frame with the corresponding taxon by readID
merge_onematch_mf <- merge(one_match_mf,speciesID_per_read_full_mf,by="readID")
# Change "ssciname" to "ID99" to reflect that this is the ID at sss>=99
colnames(merge_onematch_mf)[3] <- "ID99"

#Next, select readIDs with >1 taxonomic match and determine lowest common taxon manually
more_match_mf <- species_per_read_full_mf %>% filter(count>1)
more_match_mf$ID99 <- NA # leave this column empty for now because will fill in manually

# Manually assign lowest possible taxa to each ASV that hits to more than one taxon
# filter for each readID to see which taxa it hit 
# Also check if there is any variation in sss score
# One "tie" at genus level based on 99% sss can be further resolved by diffs in score
# Uniq5- Salvelinus alpinus and Salvelinus namaycush (assigned as S. alpinus)
speciesmf_full %>% filter(readID=="uniq22") %>% distinct(sss, ssciname)

# Add ID99s individually, based off of lowest common taxon AND
# whether or not a given taxon is known to occur in the area
more_match_mf$ID99[1] <- "Larus"
more_match_mf$ID99[2] <- "Pungitius pungitius"
more_match_mf$ID99[3] <- "Anguilla anguilla"
more_match_mf$ID99[4] <- "Felis catus"
more_match_mf$ID99[5] <- "Cottus"
more_match_mf$ID99[6] <- "Pungitius pungitius"
more_match_mf$ID99[7] <- "Alces alces"
more_match_mf$ID99[8] <- "Salvelinus fontinalis"
more_match_mf$ID99[9] <- "Cottus"
more_match_mf$ID99[10] <- "Sus scrofa"
more_match_mf$ID99[11] <- "Ursus"
more_match_mf$ID99[12] <- "Bos taurus"
more_match_mf$ID99[13] <- "Canis lupus (familiaris)"
more_match_mf$ID99[14] <- "Gadus"
more_match_mf$ID99[15] <- "Homo sapiens"
more_match_mf$ID99[16] <- "Canis lupus (familiaris)"
more_match_mf$ID99[17] <- "Esox lucius"
more_match_mf$ID99[18] <- "Pseudopleuronectes americanus"
more_match_mf$ID99[19] <- "Sus scrofa"
more_match_mf$ID99[20] <- "Stichaeidae"
more_match_mf$ID99[21] <- "Liparis gibbus"
more_match_mf$ID99[22] <- "Myoxocephalus"
more_match_mf$ID99[23] <- "Falcipennis canadensis"
more_match_mf$ID99[24] <- "Gallus gallus"
more_match_mf$ID99[25] <- "Alosa pseudoharengus"
more_match_mf$ID99[26] <- "Salvelinus alpinus"
more_match_mf$ID99[27] <- "Turdus migratorius"
more_match_mf$ID99[28] <- "Salvelinus alpinus"
more_match_mf$ID99[29] <- "Canis lupus (familiaris)"
more_match_mf$ID99[30] <- "Phoca"
more_match_mf$ID99[31] <- "Oncorhynchus gorbuscha"
more_match_mf$ID99[32] <- "Salmo trutta"

#Bind dfs with single matches and those determined manually back together
# This results in a df with a single "best" taxonomic ID for each ASV with sss>=99
fullmf <- rbind(merge_onematch_mf,more_match_mf)
# Save this for use comparing sequences later
write.csv(fullmf,"asv_species_MIFISHU.csv",row.names=FALSE)

# -------- Manipulation of ASV read counts per sample data -------------
# Load in data and format
asv_mifish_full <- read.csv(here("Compute Canada Transfers","MIFISHU","asv_counts_mifishu_full.csv"))
colnames(asv_mifish_full)[1] <- "site"
asvbysitemf_full <- as.data.frame(t(asv_mifish_full)) #transverse
colnames(asvbysitemf_full) <- asv_mifish_full$site    #rename columns with sites              
asvbysitemf_full <- asvbysitemf_full[2:671,] # get rid of first row (also site names)

asvbysitemf_full <- cbind(rownames(asvbysitemf_full), #moves row names to be first column instead
                          data.frame(asvbysitemf_full, row.names=NULL))
colnames(asvbysitemf_full)[1] <- "readID"

# Merge read counts per site by taxon ID using readID:
asvbysitemf_full <- merge(fullmf,asvbysitemf_full,by="readID", all.x=TRUE)

# Format column names so are more readable/usable (drop stuff added through bioinf pipeline):
# Results in column name just as sample names
colnames(asvbysitemf_full) <- gsub(x = colnames(asvbysitemf_full), 
                                   pattern = ".MIFISHU_merged.assembled", 
                                   replacement = "",fixed=T)
colnames(asvbysitemf_full) <- gsub(x = colnames(asvbysitemf_full), 
                                   pattern = "X12S", 
                                   replacement = "",fixed=T)
colnames(asvbysitemf_full) <- gsub(x = colnames(asvbysitemf_full), 
                                   pattern = ".", 
                                   replacement = "",fixed=T)

asvbysitemf_full <- as.data.frame(t(asvbysitemf_full)) #transverse
colnames(asvbysitemf_full) <- asvbysitemf_full[3,] #rename columns to be readIDs
asvbysitemf_full <- asvbysitemf_full[(-(1:3)),] #get rid of first 3 rows (redundant)
asvbysitemf_full<- cbind(rownames(asvbysitemf_full), 
                         data.frame(asvbysitemf_full, row.names=NULL)) #move row names to be first column
colnames(asvbysitemf_full)[1] <- "site"


# ---------Separate out LAB blank samples -----------------------
# Create separate df with just lab blanks (come back to this later/below)
blanks_full_mf <- asvbysitemf_full %>% filter(grepl("BLANK", site))
# Then separate out non-lab blanks- use this going forward
sites_full_mf <- asvbysitemf_full %>% filter(!grepl("BLANK",site))
str(sites_full_mf)
sites_full_mf[,2:101] <- sapply(sites_full_mf[,2:101],as.numeric)

#----Separate out FIELD negatives----
# Filter for neg rows
#select any row with letter N at the end of site name OR letter N followed by one number
neg_full_mf <- sites_full_mf %>%
  filter(grepl("N[0-9]$", site) | grepl("[N]$",site))  
# Create stub names for unique sites (i.e. take N off end of site name)
# Only remove last character if an N (some have numbers after, need those still)
neg_full_mf$stub <- sub("[N]$","",neg_full_mf$site)

# Select stubs that end in N#; need to trim off N but keep last number
neg_full_tocut_mf <- neg_full_mf %>% filter(grepl("N[0-9]$",stub))
# All have 5 characters with the N signifying negative being the 4th, so trim the 4th character:
neg_full_tocut_mf$stub <- sub("^(.{3}).", "\\1", neg_full_tocut_mf$stub)
# Select stubs taht do not end in N#
neg_full_nocut_mf <- neg_full_mf %>% filter(!grepl("N[0-9]$",stub))
#Bind the cut and non-cut dfs back together:
neg_full_mf <- rbind(neg_full_tocut_mf,neg_full_nocut_mf)

#copy neg_full with site names for later
neg_full_mf1 <- neg_full_mf
#Save neg_full_mf as csv
write.csv(neg_full_mf,"field_negs_MIFISHU.csv",row.names=FALSE)

#filter for non-neg rows (actual field samples)
#select rows that do not end in N AND do not end in N#
normal_full_mf <- sites_full_mf %>% filter(!grepl("N[0-9]$", site) & !grepl("[N]$",site))
# Create stub names for unique samples
# Only remove last character if a letter (some have numbers after, need those still)
normal_full_mf$stub <- sub("[A-Z]$","",normal_full_mf$site)

# For sites that have more than one negative, need to separate into corresponding samples:
# These are in neg_full_tocut_mf: HMP, LHR, MBP, PRI, RNR, SGR
normal_full_tocut_mf <- normal_full_mf %>% filter(grepl("HMP",stub) |
                                                    grepl("LHR",stub) |
                                                    grepl("MBP",stub) |
                                                    grepl("PRI",stub) |
                                                    grepl("RNR",stub) |
                                                    grepl("SGR",stub))
# All have 5 characters with the N signifying negative being the 4th, so trim the 4th character:
normal_full_tocut_mf$stub <- sub("^(.{3}).", "\\1", normal_full_tocut_mf$stub)
# Filter for ones that don't have to be trimmed this way:
normal_full_nocut_mf <- normal_full_mf %>% filter(!grepl("HMP",stub) &
                                                    !grepl("LHR",stub) &
                                                    !grepl("MBP",stub) &
                                                    !grepl("PRI",stub) &
                                                    !grepl("RNR",stub) &
                                                    !grepl("SGR",stub))

# this next part is to deal with the fact that some sample names have additional numbers at the end
# i.e. "C1", so the previous step where the only the last character was removed does not necessarily result in the same stub name as in the neg file
# Need to remove an additional letter, but only from sample names longer than 3 letters
fourplus_norm_full_mf <- normal_full_nocut_mf %>% 
  filter(nchar(stub)>3) #select sample names > 3 letters
three_norm_full_mf <- normal_full_nocut_mf %>% filter(nchar(stub) <=3) #select sample names= 3 letters
fourplus_norm_full_mf$stub <- sub("[A-Z][1-9]$","",fourplus_norm_full_mf$stub) #remove last 2 characters from sample names > 3 letters if sequence is a letter followed by number
allsites_norm_full_mf <- rbind(fourplus_norm_full_mf,three_norm_full_mf) #bind back together
allsites_norm_full_mf <- rbind(normal_full_tocut_mf,allsites_norm_full_mf)
# Previous renaming messed up OLL because the naming convention was different
allsites_norm_full_mf[216:217,102] <- "OLL"

# save as csv file
write.csv(allsites_norm_full_mf, "field_samples_MIFISHU.csv",row.names=FALSE)

#---------- Subtracting read counts in field negatives from field samples-----
#Get rid of site row (not needed anymore)
neg_full_mf$site=NULL
# Rename species columns with "neg" prefix (so can subtract from corresponding sample columns later)
names(neg_full_mf)[1:100] <- paste0("Neg_", names(neg_full_mf)[1:100])

# There are some sites for which there is no neg sample (assumed no reads?)
# Keep all samples in allsites_norm to maintain full dataset
# Convert zeros to NAs in allsites_norm df to avoid subtracting from a read count of zero (already no detection; will not be confused as a sample that is blanked out)
allsites_norm_full_mf[allsites_norm_full_mf == 0] <- NA
final_df_full_mf <- merge(allsites_norm_full_mf, neg_full_mf, 
                          by=c("stub"), all.x=TRUE)
#convert NAs to zeros for negative samples (since assume that lack of negative for site means that no reads were detected??)
final_df_full_mf[,c(103:202)][is.na(final_df_full_mf)[,c(103:202)]] <- 0

# BLOCKWISE (MULTI-COLUMN) SUBTRACTION (subtracts read count from columns with "Neg-" prefix from corresponding ASV/taxon IDs in actual samples)
final_df_full_mf[,c(3:102)] <- final_df_full_mf[,c(3:102)]-final_df_full_mf[,c(103:202)]

# REMOVE NEG SPECIES COLUMNS (dont need them anymore)
# This results in a df with field neg- corrected read counts for each sample at each site
final_df_full_mf[, grep("^Neg", names(final_df_full_mf))] <- NULL

# ****** Dealing with samples with no field blank *********
# subtract out average read depth per ASV for sites sampled on same day (at least 2 for each missing sample)
# Sites without sequenced field blanks:
# 1102019 -sampled on 29 Aug 2019; same day as NAC2019, PAL2019,1092019 (also 4 other with no field blank)
# BRG2019 -sampled on 30 Aug 2019; same day as NOR2019, BLD2019, AVA2019
# KAN2019 -sampled on 29 Aug 2019; same day as NAC2019, PAL2019,1092019 (also 4 other with no field blank)
# MKB2019 -sampled on 1 Sept 2019; same day as ENG2019, PAM2019,TOM, TOM2, TOM3
# PAN2019 -sampled on 29 Aug 2019; same day as NAC2019, PAL2019,1092019 (also 4 other with no field blank)
# SUS2019 -sampled on 27 Aug 2019; same day as NAS2019, MLR2019, DBL2019, KEN2019
# SWA2019 -sampled on 29 Aug 2019; same day as NAC2019, PAL2019,1092019 (also 4 other with no field blank)

# Start with sites sampled on 29 Aug: select from neg_full df:
neg_full_29aug_mf <- neg_full_mf %>% filter(stub=="NAC2019" | stub=="PAL2019" | 
                                              stub=="1092019")
# There are no blanks for these sites- assume blanks were clean for this marker

# Next, site sampled on 27 Aug: select from neg_full df:
neg_full_27aug_mf <- neg_full_mf %>% filter(stub=="NAS2019" | stub=="MLR2019" | 
                                              stub=="DBL2019" | stub=="KEN2019")
# take average
neg_full_27aug_mf <- data.frame(stub="mean",t(colMeans(neg_full_27aug_mf[1:100])))

# from final_df_full df, select samples with no field blank for Aug 27:
final_df_full_aug27_mf <- final_df_full_mf %>% filter(stub=="SUS2019")
# merge the two 27 Aug dfs together
aug27_mf <- cbind(final_df_full_aug27_mf[,1:2],((final_df_full_aug27_mf[,-(1:2)]) - (neg_full_27aug_mf[rep(1,times=1),-1])))

# Next, site sampled on 30 Aug: select from neg_full df:
neg_full_30aug_mf <- neg_full_mf %>% filter(stub=="NOR2019" | stub=="BLD2019" | 
                                              stub=="AVA2019")
# take average
neg_full_30aug_mf <- data.frame(stub="mean",t(colMeans(neg_full_30aug_mf[1:100])))

# from final_df_full df, select samples with no field blank for Aug 30:
final_df_full_aug30_mf <- final_df_full_mf %>% filter(stub=="BRG2019")
# merge the two 30 Aug dfs together
aug30_mf <- cbind(final_df_full_aug30_mf[,1:2],((final_df_full_aug30_mf[,-(1:2)]) - (neg_full_30aug_mf[rep(1,times=3),-1])))

# Next, site sampled on 1 Sept: select from neg_full df:
neg_full_1sept_mf <- neg_full_mf %>% filter(stub=="ENG2019" | stub=="PAM2019" | 
                                              stub=="TOM2019" | stub=="TOM22019" |
                                              stub=="TOM32019")
# take average
neg_full_1sept_mf <- data.frame(stub="mean",t(colMeans(neg_full_1sept_mf[1:100])))

# from final_df_full df, select samples with no field blank for sept1:
final_df_full_sept1_mf <- final_df_full_mf %>% filter(stub=="MKB2019")
# merge the two sept 1 dfs together
sept1_mf <- cbind(final_df_full_sept1_mf[,1:2],((final_df_full_sept1_mf[,-(1:2)]) - (neg_full_1sept_mf[rep(1,times=3),-1])))

# Now need to combine these specific sites back with full data
# Don't need to exclude Aug 29th sites because there were no field blanks for other samples that day (clean)
final_df_full_withblank_mf <- final_df_full_mf %>% filter(stub!="MKB2019" & stub!="BRG2019" &
                                                            stub!="SUS2019")
final_df_full_mf <- rbind(final_df_full_withblank_mf,aug27_mf,aug30_mf, sept1_mf)


# Save final_df_full_mf as field blank-corrected csv
write.csv(final_df_full_mf, "field_blank_corrected_MIFISHU.csv",row.names=FALSE)


#-----Subtracting reads in lab blanks from field samples---------------------
# Subtract maximum read count per ASV from each sample
blanks_full_mf[,2:101] <- sapply(blanks_full_mf[2:101],as.numeric)
blanks_full_mf1 <- rbind(blanks_full_mf, data.frame(site = 'max', 
                                                    summarise_each(blanks_full_mf[,-1], funs(max(., na.rm=TRUE)))))
final_df_full_mf[final_df_full_mf <= 0] <- NA
lab_blank_max_sub_mf <- cbind(final_df_full_mf[,1:2],((final_df_full_mf[,-(1:2)]) - (blanks_full_mf1[rep(19,times=542),-1])))

# Format data frame to be presence/absence; drop site column and rename stub as site
binarysites_full_mf_max <- lab_blank_max_sub_mf[,-2]
# filter so that singletons are discarded (i.e. presence= read count >1)
binarysites_full_mf_max2 <- as.data.frame(cbind(binarysites_full_mf_max$stub,(+(binarysites_full_mf_max[-1] > 1))))
colnames(binarysites_full_mf_max2)[1] <- "site"
colnames(binarysites_full_mf_max2) <- gsub(x = colnames(binarysites_full_mf_max2), #remove periods from column names
                                           pattern = ".", replacement = " ",fixed=T)
binarysites_full_mf_max2[,2:101] <- sapply(binarysites_full_mf_max2[2:101],as.numeric)
# save this file, will be input for number of ASVs per sample and marker
write.csv(binarysites_full_mf_max2,"mifish_asvs_persample.csv",
          row.names=FALSE)

# Aggregate samples by site
aggsite_full_mf_max <- binarysites_full_mf_max2 %>% group_by(site) %>% 
  summarise_each(funs(sum(., na.rm=TRUE))) #this sums pres/abs counts for each taxon

# a) if wanting to look at detections across field replicates:
# Manually select which columns are fish
onlyfish_mf_sep_max <- aggsite_full_mf_max[,c(1,2,3,6,8:10,15,16,21,23,
                                              26,28,30,35,36,40:42,44,46,
                                              48:50,54,55,58,59,62,63,
                                              70,74,79,80,90,91,97)]
# Gather 
library(tidyr)
onlyfish_mf_sep_max2 <- gather(onlyfish_mf_sep_max, taxon, pres_abs, 
                               c("Gasterosteus aculeatus":"Salmo trutta"), 
                               factor_key=TRUE)

# Formatting column names to combine multiple ASVs per species:
# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
onlyfish_mf_sep_max2$taxon <- sub("[0-9]$","",onlyfish_mf_sep_max2$taxon)
# Get rid of spaces at the end of taxa
onlyfish_mf_sep_max2$taxon <- sub("[ ]$","",onlyfish_mf_sep_max2$taxon)
# Merge by taxa and add pres_abs
onlyfish_mf_sep_max2$site <- as.factor(onlyfish_mf_sep_max2$site)
onlyfish_mf_sep_max2$taxon <- as.factor(onlyfish_mf_sep_max2$taxon)
onlyfish_mf_sep_max2$pres_abs <- as.numeric(onlyfish_mf_sep_max2$pres_abs)
onlyfish_mf_sep_max2<- onlyfish_mf_sep_max2 %>% group_by(site,taxon) %>% 
  summarise_each(funs(sum))
# Save as input for work classifying certainty levels based on number of field detections
write.csv(onlyfish_mf_sep_max2,"onlyfish_mf_sep_max2.csv",row.names=FALSE)

# b) If wanting to look at just pres/abs across all field reps per site:
#Convert back to just pres/abs (i.e. pres if >0, abs if 0)
binaryaggsite_full_mf_max <-as.data.frame(cbind(aggsite_full_mf_max$site,(+(aggsite_full_mf_max[-1] > 0))))
colnames(binaryaggsite_full_mf_max)[1] <- "site"
# Manually select which columns are fish
onlyfish_full_mf_max <- binaryaggsite_full_mf_max[,c(1,2,3,6,8:10,15,16,21,23,
                                                     26,28,30,35,36,40:42,44,46,
                                                     48:50,54,55,58,59,62,63,
                                                     70,74,79,80,90,91,97)]
onlyfish_full_mf_max[,2:37] <- sapply(onlyfish_full_mf_max[,2:37],as.numeric)

# Formatting column names to combine multiple ASVs per species:
library(tidyr)
# Gather 
onlyfish_full_mf_max2 <- gather(onlyfish_full_mf_max, taxon, pres_abs, 
                                c("Gasterosteus aculeatus":"Salmo trutta"), 
                                factor_key=TRUE)
# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
onlyfish_full_mf_max2$taxon <- sub("[0-9]$","",onlyfish_full_mf_max2$taxon)
# Get rid of spaces at the end of taxa
onlyfish_full_mf_max2$taxon <- sub("[ ]$","",onlyfish_full_mf_max2$taxon)
# Merge by taxa and add pres_abs
onlyfish_full_mf_max2$site <- as.factor(onlyfish_full_mf_max2$site)
onlyfish_full_mf_max2$taxon <- as.factor(onlyfish_full_mf_max2$taxon)
onlyfish_full_mf_max2$pres_abs <- as.numeric(onlyfish_full_mf_max2$pres_abs)
uniqASVfish_mf_max<- onlyfish_full_mf_max2 %>% group_by(site,taxon) %>% 
  summarise_each(funs(sum))
# Spread back out again
uniqASVfish_mf_max <- spread(uniqASVfish_mf_max, key = taxon, value = pres_abs)
# This is the point where the MIFISHU data can be combined with other marker data(???)
#save as csv
write.csv(uniqASVfish_mf_max,"uniqASVfish_maxlab_MIFISHU.csv",row.names=FALSE)


#-----Create df with read counts/taxon/field sample AND total raw read depth/field sample---------
# drop stub column
lab_blank_max_sub_mf <- lab_blank_max_sub_mf[,-1]
# filter so that singletons are discarded (i.e. presence= read count >1)
lab_blank_max_sub_mf[,2:101][lab_blank_max_sub_mf[,2:101] <2] <- 0
# Convert NAs to zeros as well
lab_blank_max_sub_mf[,2:101][is.na(lab_blank_max_sub_mf[,2:101])] <- 0

# Get total raw read depth per sample:
asv_mifish_full <- read.csv(here("Compute Canada Transfers","MIFISHU","asv_counts_mifishu_full.csv"))
colnames(asv_mifish_full)[1] <- "site"
# clean up site names
asv_mifish_full$site <- gsub("-MIFISHU_merged.assembled", "", asv_mifish_full$site)
asv_mifish_full$site <- gsub("12S-", "", asv_mifish_full$site)
asv_mifish_full$site <- gsub("-", "", asv_mifish_full$site)

asv_mifish_full$tot_depth <- rowSums(asv_mifish_full[,2:671])
#Cut to just total depth
sample_depth_mifish <- asv_mifish_full[,c(1,672)]

# Merge sample_depth df with lab_blank_sub_max df
mifish_depths <- merge(lab_blank_max_sub_mf,sample_depth_mifish,all.x=TRUE,all.y=FALSE)

# Keep only fish
mifish_fish_depths <- mifish_depths[,c(1,2,3,6,8:10,15,16,21,23,
                                       26,28,30,35,36,40:42,44,46,
                                       48:50,54,55,58,59,62,63,
                                       70,74,79,80,90,91,97,102)]
colnames(mifish_fish_depths) <- gsub(x = colnames(mifish_fish_depths), #remove periods from column names
                                     pattern = ".", replacement = " ",fixed=T)


# Formatting column names to combine multiple ASVs per species:
# Gather 
library(tidyr)
mifish_fish_depths2 <- gather(mifish_fish_depths, taxon, depth, 
                              c("Gasterosteus aculeatus":"Salmo trutta"), 
                              factor_key=TRUE)

# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
mifish_fish_depths2$taxon <- sub("[0-9]$","",mifish_fish_depths2$taxon)
# Get rid of spaces at the end of taxa
mifish_fish_depths2$taxon <- sub("[ ]$","",mifish_fish_depths2$taxon)
# Merge by taxa and add pres_abs
mifish_fish_depths2$site <- as.factor(mifish_fish_depths2$site)
mifish_fish_depths2$taxon <- as.factor(mifish_fish_depths2$taxon)
mifish_fish_depths2$depth <- as.numeric(mifish_fish_depths2$depth)
mifish_fish_depths3 <- mifish_fish_depths2 %>% group_by(site,taxon,tot_depth) %>% 
  summarise_each(funs(sum))
mifish_fish_depths3$marker <- "MIFISHU"

write.csv(mifish_fish_depths3,"mifish_depths_perfieldsample.csv",row.names=FALSE)

