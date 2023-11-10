
###########################################################################
#------Script for Processing BLAST output from FISHE marker data -----------

###########################################################################

#-------BLAST data processing------------------
# Load in BLAST output data (transferred from Compute Canada)
library(here)
blastfish_full <- read.csv(here("Compute Canada Transfers", "FISHE", "blast_FISHE_range.csv"))

# Create column for "sequence similarity score (sss) based on percent identity (pident) and percent query coverage (qcovs)
blastfish_full$sss <- (blastfish_full$pident)*(blastfish_full$qcovs)/100

# Create new data frame only selecting rows with >=99% sss
library(dplyr)
speciesfish_full <- filter(blastfish_full,sss>=99)

# Create data frame for each ASV with sss>=99 and number of different taxonomic hits ("count")
species_per_read_full_fish <- data.frame(speciesfish_full %>% group_by(readID) %>%
                                           summarise(count = n_distinct(ssciname)))

# Create data frame listing each taxonomic hit for each ASV with sss>=99
speciesID_per_read_full_fish <- data.frame(speciesfish_full %>% group_by(readID) %>%
                                             distinct(ssciname))

#----- Assign BLAST ASV IDs to reads at each site-----
# Start easy; select readIDs with only one taxonomic match
one_match_fish <- data.frame(species_per_read_full_fish %>% filter(count==1))
# Merge this data frame with the corresponding taxon by readID
merge_onematch_fish <- merge(one_match_fish,speciesID_per_read_full_fish,by="readID")
# Change "ssciname" to "ID99" to reflect that this is the ID at sss>=99
colnames(merge_onematch_fish)[3] <- "ID99"

#Next, select readIDs with >1 taxonomic match and determine lowest common taxon manually
more_match_fish <- species_per_read_full_fish %>% filter(count>1)
more_match_fish$ID99 <- NA # leave this column empty for now because will fill in manually

# Manually assign lowest possible taxa to each ASV that hits to more than one taxon
#filter for each readID to see which taxa it hit 
# Also check if there is any variation in sss score
speciesfish_full %>% filter(readID=="uniq1529") %>% 
  distinct(sss, ssciname)
# Some "ties" at genus level based on 99% sss can be further resolved by diffs in score
# Uniq1252- A. rostrata and A. anguilla- assigned as A. rostrata
# Uniq1529- C. cognatus and C. bairdii- assigned as C. cognatus
# Uniq489986- A. pseudoharengus and A. aestivalis- assigned as A. pseudoharengus
# Uniq500653- S. salar and S. trutta- assigned as S. salar
# Uniq5237- C. commersonii and C. catostomus- assigned as C. commersonii

# Add ID99s individually, based off of lowest common taxon AND
# whether or not a given taxon is known to occur in the area
more_match_fish$ID99[1] <- "Anas"
more_match_fish$ID99[2] <- "Larus"
more_match_fish$ID99[3] <- "Anaxyrus"
more_match_fish$ID99[4] <- "Anas"
more_match_fish$ID99[5] <- "Aythya"
more_match_fish$ID99[6] <- "Gallus gallus"
more_match_fish$ID99[7] <- "Anguilla rostrata"
more_match_fish$ID99[8] <- "Salvelinus alpinus"
more_match_fish$ID99[9] <- "Bos taurus"
more_match_fish$ID99[10] <- "Anaxyrus"
more_match_fish$ID99[11] <- "Larus"
more_match_fish$ID99[12] <- "Salvelinus alpinus"
more_match_fish$ID99[13] <- "Canis lupus"
more_match_fish$ID99[14] <- "Cottus cognatus"
more_match_fish$ID99[15] <- "Aves"
more_match_fish$ID99[16] <- "Cottus"
more_match_fish$ID99[17] <- "Mallotus villosus"
more_match_fish$ID99[18] <- "Myoxocephalus"
more_match_fish$ID99[19] <- "Tautogolabrus adspersus"
more_match_fish$ID99[20] <- "Salmo trutta"
more_match_fish$ID99[21] <- "Anguilla anguilla"
more_match_fish$ID99[22] <- "Boreogadus saida"
more_match_fish$ID99[23] <- "Sus scrofa"
more_match_fish$ID99[24] <- "Gadus morhua"
more_match_fish$ID99[25] <- "Anas"
more_match_fish$ID99[26] <- "Branta"
more_match_fish$ID99[27] <- "Myoxocephalus"
more_match_fish$ID99[28] <- "Bos taurus"
more_match_fish$ID99[29] <- "Phoca"
more_match_fish$ID99[30] <- "Zonotrichia leucophrys"
more_match_fish$ID99[31] <- "Homo sapiens"
more_match_fish$ID99[32] <- "Canis lupus (familiaris)"
more_match_fish$ID99[33] <- "Anatidae"
more_match_fish$ID99[34] <- "Pseudopleuronectes americanus"
more_match_fish$ID99[35] <- "Homo sapiens"
more_match_fish$ID99[36] <- "Alosa pseudoharengus"
more_match_fish$ID99[37] <- "Salmo salar"
more_match_fish$ID99[38] <- "Anguilla anguilla"
more_match_fish$ID99[39] <- "Catostomus commersonii"
more_match_fish$ID99[40] <- "Catostomus catostomus"
more_match_fish$ID99[41] <- "Rangifer tarandus"
more_match_fish$ID99[42] <- "Gasterosteus aculeatus"
more_match_fish$ID99[43] <- "Pungitius pungitius"
more_match_fish$ID99[44] <- "Felis"
more_match_fish$ID99[45] <- "Myoxocephalus"
more_match_fish$ID99[46] <- "Fundulus diaphanus"
more_match_fish$ID99[47] <- "Salvelinus fontinalis"
more_match_fish$ID99[48] <- "Lithobates clamitans"
more_match_fish$ID99[49] <- "Aythya"
more_match_fish$ID99[50] <- "Alces"


#Bind dfs with single matches and those determined manually back together
# This results in a df with a single "best" taxonomic ID for each ASV with sss>=99
fullfish <- rbind(merge_onematch_fish,more_match_fish)
# Save this for use comparing sequences later
write.csv(fullfish,"asv_species_FISHE.csv",row.names=FALSE)

# -------- Manipulation of ASV read counts per sample data -------------
# Load in data and format
asv_fish_full <- read.csv(here("Compute Canada Transfers","FISHE",
                               "reduced_asv_counts.csv"))

# Merge read counts per site by taxon ID using readID:
asvbysite_full_fish <- merge(fullfish,asv_fish_full,by=c("readID","count"), 
                             all.x=TRUE)

# Get rid of count column
asvbysite_full_fish <- asvbysite_full_fish[,-2]

# Format column names so are more readable/usable (drop stuff added through bioinf pipeline):
# Results in column name just as sample names
colnames(asvbysite_full_fish) <- gsub(x = colnames(asvbysite_full_fish), 
                                      pattern = "FISHE_merged.assembled", replacement = "",fixed=T)
colnames(asvbysite_full_fish) <- gsub(x = colnames(asvbysite_full_fish), 
                                      pattern = "X12S", replacement = "",fixed=T)
colnames(asvbysite_full_fish) <- gsub(x = colnames(asvbysite_full_fish), 
                                      pattern = ".", replacement = "",fixed=T)

asvbysite_full_fish <- as.data.frame(t(asvbysite_full_fish)) #transverse
colnames(asvbysite_full_fish) <- asvbysite_full_fish[2,] #rename columns to be readIDs
asvbysite_full_fish <- asvbysite_full_fish[(-(1:2)),] #get rid of first 3 rows (redundant)
asvbysite_full_fish<- cbind(rownames(asvbysite_full_fish), 
                            data.frame(asvbysite_full_fish, row.names=NULL)) #move row names to be first column
colnames(asvbysite_full_fish)[1] <- "site"

# ---------Separate out LAB blank samples -----------------------
# Create separate df with just lab blanks (come back to this later/below)
blanks_full_fish <- asvbysite_full_fish %>% filter(grepl("BLANK", site))
# Then separate out non-lab blanks- use this going forward
sites_full_fish <- asvbysite_full_fish %>% filter(!grepl("BLANK",site))
str(sites_full_fish)
sites_full_fish[,2:121] <- sapply(sites_full_fish[,2:121],as.numeric)

#----Separate out FIELD negatives----
# Filter for neg rows
#select any row with letter N at the end of site name OR letter N followed by one number
neg_full_fish <- sites_full_fish %>% 
  filter(grepl("N[0-9]$", site) | grepl("[N]$",site)) 
# Create stub names for unique samples (i.e. take N off end of site name)
# Only remove last character if an N (some have numbers after, need those still)
neg_full_fish$stub <- sub("[N]$","",neg_full_fish$site)

# Select stubs that end in N#; need to trim off N but keep last number
neg_full_tocut_fish <- neg_full_fish %>% filter(grepl("N[0-9]$",stub))
# All have 5 characters with the N signifying negative being the 4th, so trim the 4th character:
neg_full_tocut_fish$stub <- sub("^(.{3}).", "\\1", neg_full_tocut_fish$stub)
# Select stubs taht do not end in N#
neg_full_nocut_fish <- neg_full_fish %>% filter(!grepl("N[0-9]$",stub))
#Bind the cut and non-cut dfs back together:
neg_full_fish <- rbind(neg_full_tocut_fish,neg_full_nocut_fish)

# copy neg_full_fish for later
neg_full_fish1 <- neg_full_fish
#Save neg_full as csv
write.csv(neg_full_fish,"field_negs_FISHE.csv",row.names=FALSE)

#filter for non-neg rows (actual field samples)
#select rows that do not end in N AND do not end in N#
normal_full_fish <- sites_full_fish %>% 
  filter(!grepl("N[0-9]$", site) & !grepl("[N]$",site))
# Create stub names for unique samples
# Only remove last character if a letter (some have numbers after, need those still)
normal_full_fish$stub <- sub("[A-Z]$","",normal_full_fish$site)

# For sites that have more than one negative, need to separate into corresponding samples:
# These are in neg_full_tocut_fish: BVR, CBC, HDB, HMP, LHR, MBP, PRI, PRR, RNR,SGR
normal_full_tocut_fish <- normal_full_fish %>% filter(grepl("BVR",stub) |
                                                        grepl("CBC",stub) |
                                                        grepl("HDB",stub) |
                                                        grepl("HMP",stub) |
                                                        grepl("LHR",stub) |
                                                        grepl("MBP",stub) |
                                                        grepl("PRI",stub) |
                                                        grepl("PRR",stub) |
                                                        grepl("RNR",stub) |
                                                        grepl("SGR",stub))
# All have 5 characters with the N signifying negative being the 4th, so trim the 4th character:
normal_full_tocut_fish$stub <- sub("^(.{3}).", "\\1", normal_full_tocut_fish$stub)
# Filter for ones that don't have to be trimmed this way:
normal_full_nocut_fish <- normal_full_fish %>% filter(!grepl("BVR",stub) &
                                                        !grepl("CBC",stub) &
                                                        !grepl("HDB",stub) &
                                                        !grepl("HMP",stub) &
                                                        !grepl("LHR",stub) &
                                                        !grepl("MBP",stub) &
                                                        !grepl("PRI",stub) &
                                                        !grepl("PRR",stub) &
                                                        !grepl("RNR",stub) &
                                                        !grepl("SGR",stub))

# this next part is to deal with the fact that some sample names have additional numbers at the end
# i.e. "C1", so the previous step where the only the last character was removed does not necessarily result in the same stub name as in the neg file
# Need to remove an additional letter, but only from sample names longer than 3 letters
fourplus_norm_full_fish <- normal_full_nocut_fish %>% filter(nchar(stub)>3) #select sample names > 3 letters
three_norm_full_fish <- normal_full_nocut_fish %>% filter(nchar(stub) <=3) #select sample names= 3 letters
fourplus_norm_full_fish$stub <- sub("[A-Z][1-9]$","",fourplus_norm_full_fish$stub) #remove last 2 characters from sample names > 3 letters if sequence is a letter followed by number
allsites_norm_full_fish <- rbind(fourplus_norm_full_fish,three_norm_full_fish) #bind back together
allsites_norm_full_fish <- rbind(normal_full_tocut_fish,allsites_norm_full_fish)
# Previous renaming messed up OLL because the naming convention was different
allsites_norm_full_fish[227:229,122] <- "OLL"

#save as CSV
write.csv(allsites_norm_full_fish,"field_samples_FISHE.csv",row.names=FALSE)

# --------Subtracting read counts in field negatives from field samples-------------
#Get rid of site row (not needed anymore)
neg_full_fish$site=NULL
# Rename species columns with "neg" prefix (so can subtract from corresponding sample columns later)
names(neg_full_fish)[1:120] <- paste0("Neg_", names(neg_full_fish)[1:120])

# There are some sites for which there is no neg sample (assumed no reads?)
# Keep all samples in allsites_norm to maintain full dataset
# Convert zeros to NAs in allsites_norm df to avoid subtracting from a read count of zero (already no detection; will not be confused as a sample that is blanked out)
allsites_norm_full_fish[allsites_norm_full_fish == 0] <- NA
final_df_full_fish <- merge(allsites_norm_full_fish, neg_full_fish, 
                            by=c("stub"), all.x=TRUE)
#convert NAs to zeros for negative samples (since assume that lack of negative for site means that no reads were detected??)
final_df_full_fish[,c(123:242)][is.na(final_df_full_fish)[,c(123:242)]] <- 0

# BLOCKWISE (MULTI-COLUMN) SUBTRACTION (subtracts read count from columns with "Neg-" prefix from corresponding ASV/taxon IDs in actual samples)
final_df_full_fish[,c(3:122)] <- final_df_full_fish[,c(3:122)]-final_df_full_fish[,c(123:242)]

# REMOVE NEG SPECIES COLUMNS (dont need them anymore)
# This results in a df with field neg- corrected read counts for each sample at each site
final_df_full_fish[, grep("^Neg", names(final_df_full_fish))] <- NULL

# ****** Dealing with samples with no field blank*****
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
neg_full_29aug_fish <- neg_full_fish %>% filter(stub=="NAC2019" | stub=="PAL2019" | 
                                                  stub=="1092019")
# take average
neg_full_29aug_fish <- data.frame(stub="mean",t(colMeans(neg_full_29aug_fish[1:120])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_aug29_fish <- final_df_full_fish %>% filter(stub=="1102019" | stub=="KAN2019" |
                                                            stub=="PAN2019" | stub=="SWA2019")
# merge the two 29 Aug dfs together
aug29_fish <- cbind(final_df_full_aug29_fish[,1:2],((final_df_full_aug29_fish[,-(1:2)]) - (neg_full_29aug_fish[rep(1,times=11),-1])))

# Next, site sampled on 27 Aug: select from neg_full df:
neg_full_27aug_fish <- neg_full_fish %>% filter(stub=="NAS2019" | stub=="MLR2019" | 
                                                  stub=="DBL2019" | stub=="KEN2019")
# take average
neg_full_27aug_fish <- data.frame(stub="mean",t(colMeans(neg_full_27aug_fish[1:120])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_aug27_fish <- final_df_full_fish %>% filter(stub=="SUS2019")
# merge the two 27 Aug dfs together
aug27_fish <- cbind(final_df_full_aug27_fish[,1:2],((final_df_full_aug27_fish[,-(1:2)]) - (neg_full_27aug_fish[rep(1,times=3),-1])))

# Next, site sampled on 30 Aug: select from neg_full df:
neg_full_30aug_fish <- neg_full_fish %>% filter(stub=="NOR2019" | stub=="BLD2019" | 
                                                  stub=="AVA2019")
# take average
neg_full_30aug_fish <- data.frame(stub="mean",t(colMeans(neg_full_30aug_fish[1:120])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_aug30_fish <- final_df_full_fish %>% filter(stub=="BRG2019")
# merge the two 30 Aug dfs together
aug30_fish <- cbind(final_df_full_aug30_fish[,1:2],((final_df_full_aug30_fish[,-(1:2)]) - (neg_full_30aug_fish[rep(1,times=3),-1])))

# Next, site sampled on 1 Sept: select from neg_full df:
neg_full_1sept_fish <- neg_full_fish %>% filter(stub=="ENG2019" | stub=="PAM2019" | 
                                                  stub=="TOM2019" | stub=="TOM22019" |
                                                  stub=="TOM32019")
# take average
neg_full_1sept_fish <- data.frame(stub="mean",t(colMeans(neg_full_1sept_fish[1:120])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_sept1_fish <- final_df_full_fish %>% filter(stub=="MKB2019")
# merge the two sept 1 dfs together
sept1_fish <- cbind(final_df_full_sept1_fish[,1:2],((final_df_full_sept1_fish[,-(1:2)]) - (neg_full_1sept_fish[rep(1,times=3),-1])))

# Now need to combine these specific sites back with full data
final_df_full_withblank_fish <- final_df_full_fish %>% filter(stub!="MKB2019" & stub!="BRG2019" &
                                                                stub!="SUS2019" & stub!="1102019" &
                                                                stub!="KAN2019" & stub!="PAN2019" &
                                                                stub!="SWA2019")
final_df_full_fish <- rbind(final_df_full_withblank_fish, aug27_fish,
                            aug29_fish,aug30_fish, sept1_fish)



# Save final_df_full_fish as field blank-corrected csv
write.csv(final_df_full_fish,"field_blank_corrected_FISHE.csv",row.names=FALSE)


#-----Subtracting reads in lab blanks from field samples---------------------
# Subtract maximum read count per ASV from each sample
blanks_full_fish[,2:121] <- sapply(blanks_full_fish[2:121],as.numeric)
blanks_full_fish1 <- rbind(blanks_full_fish, data.frame(site = 'max', 
                                                        summarise_each(blanks_full_fish[,-1], funs(max(., na.rm=TRUE)))))

lab_blank_sub_fish_max <- cbind(final_df_full_fish[,1:2],((final_df_full_fish[,-(1:2)]) - (blanks_full_fish1[rep(29,times=568),-1])))

# Format data frame to be presence/absence; drop site column and rename stub as site
binarysites_full_fish_max <- lab_blank_sub_fish_max[,-2]
# filter so that singletons are discarded (i.e. presence= read count >1)
binarysites_full_fish_max2 <- as.data.frame(cbind(binarysites_full_fish_max$stub,(+(binarysites_full_fish_max[-1] > 1))))
colnames(binarysites_full_fish_max2)[1] <- "site"
colnames(binarysites_full_fish_max2) <- gsub(x = colnames(binarysites_full_fish_max2), #remove periods from column names
                                             pattern = ".", replacement = " ",fixed=T)
binarysites_full_fish_max2[,2:121] <- sapply(binarysites_full_fish_max2[2:121],as.numeric)
# save this file, will be input for number of asvs per sample and marker
write.csv(binarysites_full_fish_max2,"fishe_asvs_persample.csv",
          row.names=FALSE)

# Aggregate samples by site
aggsite_full_fish_max <- binarysites_full_fish_max2 %>% group_by(site) %>% 
  summarise_each(funs(sum(., na.rm=TRUE))) #this sums pres/abs counts for each taxon

# a) if wanting to look at detections across field replicates:
# Manually select which columns are fish
onlyfish_fish_sep_max <- aggsite_full_fish_max[,c(1,12,13,17,21,24,26,28:31,
                                                  33,35,36,38,39,42,44,45,50:52,
                                                  55:57,61,73,75,78,80,82,84,
                                                  87:90,92,94,98,100,103,
                                                  106,110,111,121)]

# Gather 
library(tidyr)
onlyfish_fish_sep_max2 <- gather(onlyfish_fish_sep_max, taxon, pres_abs, 
                                 c("Anguilla rostrata":"Esox lucius 1"), 
                                 factor_key=TRUE)

# Formatting column names to combine multiple ASVs per species:
# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
onlyfish_fish_sep_max2$taxon <- sub("[0-9]$","",onlyfish_fish_sep_max2$taxon)
# Get rid of spaces at the end of taxa
onlyfish_fish_sep_max2$taxon <- sub("[ ]$","",onlyfish_fish_sep_max2$taxon)
# Merge by taxa and add pres_abs
onlyfish_fish_sep_max2$site <- as.factor(onlyfish_fish_sep_max2$site)
onlyfish_fish_sep_max2$taxon <- as.factor(onlyfish_fish_sep_max2$taxon)
onlyfish_fish_sep_max2$pres_abs <- as.numeric(onlyfish_fish_sep_max2$pres_abs)
onlyfish_fish_sep_max2<- onlyfish_fish_sep_max2 %>% group_by(site,taxon) %>% 
  summarise_each(funs(sum))
# Save as input for work classifying certainty levels based on number of field detections
write.csv(onlyfish_fish_sep_max2,"onlyfish_fish_sep_max2.csv",row.names=FALSE)


# b) If wanting to look at just pres/abs across all field reps per site:
#Convert back to just pres/abs (i.e. pres if >0, abs if 0)
binaryaggsite_full_fish_max <-as.data.frame(cbind(aggsite_full_fish_max$site,(+(aggsite_full_fish_max[-1] > 0))))
colnames(binaryaggsite_full_fish_max)[1] <- "site"
# Manually select which columns are fish
onlyfish_full_fish_max <- binaryaggsite_full_fish_max[,c(1,12,13,17,21,24,26,28:31,
                                                         33,35,36,38,39,42,44,45,50:52,
                                                         55:57,61,73,75,78,80,82,84,
                                                         87:90,92,94,98,100,103,
                                                         106,110,111,121)]

# Formatting column names to combine multiple ASVs per species:
library(tidyr)
# Gather 
onlyfish_full_fish_max2 <- gather(onlyfish_full_fish_max, taxon, pres_abs, 
                                  c("Anguilla rostrata":"Esox lucius 1"), 
                                  factor_key=TRUE)
# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
onlyfish_full_fish_max2$taxon <- sub("[0-9]$","",onlyfish_full_fish_max2$taxon)
# Get rid of spaces at the end of taxa
onlyfish_full_fish_max2$taxon <- sub("[ ]$","",onlyfish_full_fish_max2$taxon)
# Merge by taxa and add pres_abs
onlyfish_full_fish_max2$site <- as.factor(onlyfish_full_fish_max2$site)
onlyfish_full_fish_max2$taxon <- as.factor(onlyfish_full_fish_max2$taxon)
onlyfish_full_fish_max2$pres_abs <- as.numeric(onlyfish_full_fish_max2$pres_abs)
uniqASVfish_fish_max<- onlyfish_full_fish_max2 %>% group_by(site,taxon) %>% 
  summarise_each(funs(sum))
# Spread back out again
uniqASVfish_fish_max <- spread(uniqASVfish_fish_max, key = taxon, value = pres_abs)
# This is the point where the 12Steleo data can be combined with other marker data(???)
write.csv(uniqASVfish_fish_max,"uniqASVfish_maxlab_FISHE.csv",row.names=FALSE)


#-----Create df with read counts/taxon/field sample AND total raw read depth/field sample---------
# drop stub column
lab_blank_sub_fish_max <- lab_blank_sub_fish_max[,-1]
# filter so that singletons are discarded (i.e. presence= read count >1)
lab_blank_sub_fish_max[,2:121][lab_blank_sub_fish_max[,2:121] <2] <- 0
# Convert NAs to zeros as well
lab_blank_sub_fish_max[,2:121][is.na(lab_blank_sub_fish_max[,2:121])] <- 0

# Get total raw read depth per sample:
asv_fish_full <- read.csv(here("Compute Canada Transfers","FISHE",
                               "reduced_asv_counts.csv"))
asv_fish_full1 <- asv_fish_full[,-2]
asv_fish_full1 <- as.data.frame(t(asv_fish_full1)) #transverse
asv_fish_full1 <- cbind(rownames(asv_fish_full1), #moves row names to be first column instead
                        data.frame(asv_fish_full1, row.names=NULL))
colnames(asv_fish_full1) <- asv_fish_full1[1,]    #rename columns with sites              
asv_fish_full1 <- asv_fish_full1[2:752,] # get rid of first row (also site names)

colnames(asv_fish_full1)[1] <- "site"
# clean up site names
asv_fish_full1$site <- gsub(".FISHE_merged.assembled", "", asv_fish_full1$site)
asv_fish_full1$site <- gsub("X12S.", "", asv_fish_full1$site)
asv_fish_full1$site <- gsub("\\.", "", asv_fish_full1$site)

# Convert read depths to numeric and rum depths across rows (sites):
asv_fish_full1[,2:121] <- sapply(asv_fish_full1[,2:121],as.numeric)
asv_fish_full1$tot_depth <- rowSums(asv_fish_full1[,2:121])
#Cut to just total depth
sample_depth_fishe <- asv_fish_full1[,c(1,122)]

# Merge sample_depth df with lab_blank_sub_max df
fishe_depths <- merge(lab_blank_sub_fish_max,sample_depth_fishe,all.x=TRUE,all.y=FALSE)

# Keep only fish
fishe_fish_depths <- fishe_depths[,c(1,12,13,17,21,24,26,28:31,
                                     33,35,36,38,39,42,44,45,50:52,
                                     55:57,61,73,75,78,80,82,84,
                                     87:90,92,94,98,100,103,
                                     106,110,111,121,122)]

colnames(fishe_fish_depths) <- gsub(x = colnames(fishe_fish_depths), #remove periods from column names
                                    pattern = ".", replacement = " ",fixed=T)


# Formatting column names to combine multiple ASVs per species:
# Gather 
library(tidyr)
fishe_fish_depths2 <- gather(fishe_fish_depths, taxon, depth, 
                             c("Anguilla rostrata":"Esox lucius 1"), 
                             factor_key=TRUE)

# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
fishe_fish_depths2$taxon <- sub("[0-9]$","",fishe_fish_depths2$taxon)
# Get rid of spaces at the end of taxa
fishe_fish_depths2$taxon <- sub("[ ]$","",fishe_fish_depths2$taxon)
# Merge by taxa and add pres_abs
fishe_fish_depths2$site <- as.factor(fishe_fish_depths2$site)
fishe_fish_depths2$taxon <- as.factor(fishe_fish_depths2$taxon)
fishe_fish_depths2$depth <- as.numeric(fishe_fish_depths2$depth)
fishe_fish_depths3 <- fishe_fish_depths2 %>% group_by(site,taxon,tot_depth) %>% 
  summarise_each(funs(sum))
fishe_fish_depths3$marker <- "FISHE"

write.csv(fishe_fish_depths3,"fishe_depths_perfieldsample.csv",row.names=FALSE)



