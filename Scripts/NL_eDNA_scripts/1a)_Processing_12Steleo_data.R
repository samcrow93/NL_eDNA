
###########################################################################
#------Script for Processing BLAST output from 12Steleo marker data -----------

###########################################################################

# ---------BLAST data processing---------------------
# Load in BLAST output data (transferred from Compute Canada)
library(here)
blast12s_full <- read.csv(here("Compute Canada Transfers", "12Steleo", "12steleo_blast_full.csv"))

# Create column for "sequence similarity score (sss) based on percent identity (pident) and percent query coverage (qcovs)
blast12s_full$sss <- (blast12s_full$pident)*(blast12s_full$qcovs)/100

# Create new data frame only selecting rows with >=99% sss
library(dplyr)
species12s_full <- filter(blast12s_full,sss>=99)

# Create data frame for each ASV with sss>=99 and number of different taxonomic hits ("count")
species_per_read_full <- data.frame(species12s_full %>% group_by(readID) %>%
  summarise(count = n_distinct(ssciname)))

# Create data frame listing each taxonomic hit for each ASV with sss>=99
speciesID_per_read_full <- data.frame(species12s_full %>% group_by(readID) %>%
  distinct(ssciname))

#----- Assign BLAST ASV IDs to reads at each site-----
# Start easy; select readIDs with only one taxonomic match
one_match <- data.frame(species_per_read_full %>% filter(count==1))
# Merge this data frame with the corresponding taxon by readID
merge_onematch <- merge(one_match,speciesID_per_read_full,by="readID")
# Change "ssciname" to "ID99" to reflect that this is the ID at sss>=99
colnames(merge_onematch)[3] <- "ID99"

#Next, select readIDs with >1 taxonomic match and determine lowest common taxon manually
more_match <- species_per_read_full %>% filter(count>1)
more_match$ID99 <- NA # leave this column empty for now because will fill in manually

# Manually assign lowest possible taxa to each ASV that hits to more than one taxon
# do this by filtering for each readID to see which taxa it hit 
# Also check if there is any variation in sss score
species12s_full %>% filter(readID=="uniq1") %>% distinct(sss, ssciname) # run this for each relevant readID
# all of the ASVs with multiple hits have an sss score of 100
# therefore no ability to sort/pick best based on score

# Add ID99s individually, based off of lowest common taxon AND
# whether or not a given taxon is known to occur in the area
more_match$ID99[1] <- "Salvelinus"
more_match$ID99[2] <- "Anguilla rostrata"
more_match$ID99[3] <- "Gadidae"
more_match$ID99[4] <- "Lithobates clamitans"
more_match$ID99[5] <- "Catostomus commersonii"
more_match$ID99[6] <- "Gasterosteus aculeatus"
more_match$ID99[7] <- "Larus"
more_match$ID99[8] <- "Lotidae"
more_match$ID99[9] <- "Stichaeidae"
more_match$ID99[10] <- "Eumesogrammus praecisus"
more_match$ID99[11] <- "Ammodytes"
more_match$ID99[12] <- "Pungitius pungitius"
more_match$ID99[13] <- "Anas crecca"
more_match$ID99[14] <- "Myoxocephalus"
more_match$ID99[15] <- "Anas platyrhynchos"
more_match$ID99[16] <- "Anaxyrus"
more_match$ID99[17] <- "Gallus gallus"
more_match$ID99[18] <- "Cottus"
more_match$ID99[19] <- "Salmo trutta"

#Bind dfs with single matches and those determined manually back together
# This results in a df with a single "best" taxonomic ID for each ASV with sss>=99
full12s <- rbind(merge_onematch,more_match)
# Save this for use comparing sequences later
write.csv(full12s,"asv_species_12steleo.csv",row.names=FALSE)

# -------- Manipulation of ASV read counts per sample data -------------
# Load in data and format
asv_12s_full <- read.csv(here("Compute Canada Transfers","12Steleo","asv_counts_12steleo_full.csv"))
colnames(asv_12s_full)[1] <- "site"
asvbysite_full <- as.data.frame(t(asv_12s_full)) #transverse
colnames(asvbysite_full) <- asv_12s_full$site    #rename columns with sites              
asvbysite_full <- asvbysite_full[2:7869,] # get rid of first row (also site names)

asvbysite_full <- cbind(rownames(asvbysite_full), #moves row names to be first column instead
                        data.frame(asvbysite_full, row.names=NULL))
colnames(asvbysite_full)[1] <- "readID"

# Merge read counts per site by taxon ID using readID:
asvbysite_full <- merge(full12s,asvbysite_full,by="readID", all.x=TRUE)

# Format column names so are more readable/usable (drop stuff added through bioinf pipeline):
# Results in column name just as sample names
colnames(asvbysite_full) <- gsub(x = colnames(asvbysite_full), 
                        pattern = "12Steleo_merged.assembled", replacement = "",fixed=T)
colnames(asvbysite_full) <- gsub(x = colnames(asvbysite_full), 
                                 pattern = "X12S", replacement = "",fixed=T)
colnames(asvbysite_full) <- gsub(x = colnames(asvbysite_full), 
                                 pattern = ".", replacement = "",fixed=T)

asvbysite_full <- as.data.frame(t(asvbysite_full)) #transverse
colnames(asvbysite_full) <- asvbysite_full[3,] #rename columns to be readIDs
asvbysite_full <- asvbysite_full[(-(1:3)),] #get rid of first 3 rows (redundant)
asvbysite_full<- cbind(rownames(asvbysite_full), 
                       data.frame(asvbysite_full, row.names=NULL)) #move row names to be first column
colnames(asvbysite_full)[1] <- "site"

# ---------Separate out LAB blank samples -----------------------
# Create separate df with just lab blanks (come back to this later/below)
blanks_full <- asvbysite_full %>% filter(grepl("BLANK", site))
# Then separate out non-lab blanks- use this going forward
sites_full <- asvbysite_full %>% filter(!grepl("BLANK",site))
str(sites_full)
sites_full[,2:57] <- sapply(sites_full[,2:57],as.numeric)

#----Separate out FIELD negatives----
# Filter for neg rows
#select any row with letter N at the end of site name OR letter N followed by one number
neg_full <- sites_full %>% 
  filter(grepl("N[0-9]$", site) | grepl("[N]$",site)) 
# Create stub names for unique samples (i.e. take N off end of site name)
# Only remove last character if an N (some have numbers after, need those still)
neg_full$stub <- sub("[N]$","",neg_full$site)
 
# Select stubs that end in N#; need to trim off N but keep last number
neg_full_tocut <- neg_full %>% filter(grepl("N[0-9]$",stub))
# All have 5 characters with the N signifying negative being the 4th, so trim the 4th character:
neg_full_tocut$stub <- sub("^(.{3}).", "\\1", neg_full_tocut$stub)
# Select stubs taht do not end in N#
neg_full_nocut <- neg_full %>% filter(!grepl("N[0-9]$",stub))
#Bind the cut and non-cut dfs back together:
neg_full <- rbind(neg_full_tocut,neg_full_nocut)

#copy neg_full with site names for later
neg_full1 <- neg_full
#Save neg_full as csv
write.csv(neg_full,"field_negs_12steleo.csv",row.names=FALSE)

#filter for non-neg rows (actual field samples)
#select rows that do not end in N AND do not end in N#
normal_full <- sites_full %>% filter(!grepl("N[0-9]$", site) & !grepl("[N]$",site))
# Create stub names for unique samples
# Only remove last character if a letter (some have numbers after, need those still)
normal_full$stub <- sub("[A-Z]$","",normal_full$site)

# For sites that have more than one negative, need to separate into corresponding samples:
# These are in neg_full_tocut: BVR, CBC, HDB, HMP, LHR, PRI, PRR, RNR
normal_full_tocut <- normal_full %>% filter(grepl("BVR",stub) |
                                            grepl("CBC",stub) |
                                            grepl("HDB",stub) |
                                            grepl("HMP",stub) |
                                            grepl("LHR",stub) |
                                            grepl("PRI",stub) |
                                            grepl("PRR",stub) |
                                            grepl("RNR",stub))
# All have 5 characters with the N signifying negative being the 4th, so trim the 4th character:
normal_full_tocut$stub <- sub("^(.{3}).", "\\1", normal_full_tocut$stub)
# Filter for ones that don't have to be trimmed this way:
normal_full_nocut <- normal_full %>% filter(!grepl("BVR",stub) &
                                              !grepl("CBC",stub) &
                                              !grepl("HDB",stub) &
                                              !grepl("HMP",stub) &
                                              !grepl("LHR",stub) &
                                              !grepl("PRI",stub) &
                                              !grepl("PRR",stub) &
                                              !grepl("RNR",stub))

# this next part is to deal with the fact that some sample names have additional numbers at the end
# i.e. "C1", so the previous step where the only the last character was removed does not necessarily result in the same stub name as in the neg file
# Need to remove an additional letter, but only from sample names longer than 3 letters
fourplus_norm_full <- normal_full_nocut %>% filter(nchar(stub)>3) #select sample names > 3 letters
three_norm_full <- normal_full_nocut %>% filter(nchar(stub) <=3) #select sample names= 3 letters
fourplus_norm_full$stub <- sub("[A-Z][1-9]$","",fourplus_norm_full$stub) #remove last 2 characters from sample names > 3 letters if sequence is a letter followed by number
allsites_norm_full <- rbind(fourplus_norm_full,three_norm_full) #bind back together
allsites_norm_full <- rbind(normal_full_tocut,allsites_norm_full)
# Previous renaming messed up OLL because the naming convention was different
allsites_norm_full[234:236,58] <- "OLL"

#save as csv file
write.csv(allsites_norm_full,"field_samples_12steleo.csv",row.names=FALSE)

#-----Subtracting read counts in field negatives from field samples----
#Get rid of site row in neg_full file (not needed anymore)
neg_full$site=NULL
# Rename species columns with "neg" prefix (so can subtract from corresponding sample columns later)
names(neg_full)[1:56] <- paste0("Neg_", names(neg_full)[1:56])

# There are some sites for which there is no neg sample (assumed no reads/was actually blank)
# Keep all samples in allsites_norm to maintain full dataset
# Convert zeros to NAs in allsites_norm df to avoid subtracting from a read count of zero (already no detection; will not be confused as a sample that is blanked out)
allsites_norm_full[allsites_norm_full == 0] <- NA
final_df_full <- merge(allsites_norm_full, neg_full, 
                  by=c("stub"), all.x=TRUE)
#convert NAs to zeros for negative samples (since assume that lack of negative for site means negative was actually blank)
final_df_full[,c(59:114)][is.na(final_df_full)[,c(59:114)]] <- 0

# BLOCKWISE (MULTI-COLUMN) SUBTRACTION (subtracts read count from columns with "Neg-" prefix from corresponding ASV/taxon IDs in actual samples)
final_df_full[,c(3:58)] <- final_df_full[,c(3:58)]-final_df_full[,c(59:114)]

# REMOVE NEG SPECIES COLUMNS (dont need them anymore)
# This results in a df with field neg- corrected read counts for each sample at each site
final_df_full[, grep("^Neg", names(final_df_full))] <- NULL

# ******Dealing with samples with no field blank***********
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
neg_full_29aug <- neg_full %>% filter(stub=="NAC2019" | stub=="PAL2019" | 
                                        stub=="1092019")
# take average
neg_full_29aug <- data.frame(stub="mean",t(colMeans(neg_full_29aug[1:56])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_aug29 <- final_df_full %>% filter(stub=="1102019" | stub=="KAN2019" |
                                                  stub=="PAN2019" | stub=="SWA2019")
# merge the two 29 Aug dfs together
aug29 <- cbind(final_df_full_aug29[,1:2],((final_df_full_aug29[,-(1:2)]) - (neg_full_29aug[rep(1,times=12),-1])))

# Next, site sampled on 27 Aug: select from neg_full df:
neg_full_27aug <- neg_full %>% filter(stub=="NAS2019" | stub=="MLR2019" | 
                                        stub=="DBL2019" | stub=="KEN2019")
# take average
neg_full_27aug <- data.frame(stub="mean",t(colMeans(neg_full_27aug[1:56])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_aug27 <- final_df_full %>% filter(stub=="SUS2019")
# merge the two 27 Aug dfs together
aug27 <- cbind(final_df_full_aug27[,1:2],((final_df_full_aug27[,-(1:2)]) - (neg_full_27aug[rep(1,times=3),-1])))

# Next, site sampled on 30 Aug: select from neg_full df:
neg_full_30aug <- neg_full %>% filter(stub=="NOR2019" | stub=="BLD2019" | 
                                        stub=="AVA2019")
# take average
neg_full_30aug <- data.frame(stub="mean",t(colMeans(neg_full_30aug[1:56])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_aug30 <- final_df_full %>% filter(stub=="BRG2019")
# merge the two 30 Aug dfs together
aug30 <- cbind(final_df_full_aug30[,1:2],((final_df_full_aug30[,-(1:2)]) - (neg_full_30aug[rep(1,times=3),-1])))

# Next, site sampled on 1 Sept: select from neg_full df:
neg_full_1sept <- neg_full %>% filter(stub=="ENG2019" | stub=="PAM2019" | 
                                        stub=="TOM2019" | stub=="TOM22019" |
                                        stub=="TOM32019")
# take average
neg_full_1sept <- data.frame(stub="mean",t(colMeans(neg_full_1sept[1:56])))

# from final_df_full df, select samples with no field blank for Aug 29:
final_df_full_sept1 <- final_df_full %>% filter(stub=="MKB2019")
# merge the two sept 1 dfs together
sept1 <- cbind(final_df_full_sept1[,1:2],((final_df_full_sept1[,-(1:2)]) - (neg_full_1sept[rep(1,times=3),-1])))

# Now need to combine these specific sites back with full data
final_df_full_withblank <- final_df_full %>% filter(stub!="MKB2019" & stub!="BRG2019" &
                                                      stub!="SUS2019" & stub!="1102019" &
                                                      stub!="KAN2019" & stub!="PAN2019" &
                                                      stub!="SWA2019")
final_df_full <- rbind(final_df_full_withblank, aug27,aug29,aug30, sept1)


# Save final_df_full as field blank-corrected csv
write.csv(final_df_full,"field_blank_corrected_12steleo.csv",row.names=FALSE)


#----- Subtracting reads in lab blanks from field samples---------------------
# Subtract maximum read count per ASV from each sample
blanks_full[,2:57] <- sapply(blanks_full[2:57],as.numeric)
blanks_full1 <- rbind(blanks_full, data.frame(site = 'mean', 
                                              t(colMeans(blanks_full[-1]))))
blanks_full2 <- rbind(blanks_full1, data.frame(site = 'max', 
        summarise_each(blanks_full1[,-1], funs(max(., na.rm=TRUE)))))

# convert anything < 0 in final_df_full to NA (i.e. no detection after subtraction of field neg reads)
# Then will be able to look for any new negative values as ASVs that were blanked out by lab sample subtraction
final_df_full1 <- final_df_full
final_df_full1[final_df_full1 <= 0] <- NA
lab_blank_max_sub <- cbind(final_df_full1[,1:2],((final_df_full1[,-(1:2)]) - (blanks_full2[rep(24,times=595),-1])))

# Format data frame to be presence/absence; drop site column and rename stub as site
binarysites_full_max <- lab_blank_max_sub[,-2]
# filter so that singletons are discarded (i.e. presence= read count >1)
binarysites_full_max2 <- as.data.frame(cbind(binarysites_full_max$stub,(+(binarysites_full_max[-1] > 1))))
colnames(binarysites_full_max2)[1] <- "site"
colnames(binarysites_full_max2) <- gsub(x = colnames(binarysites_full_max2), #remove periods from column names
                                    pattern = ".", replacement = " ",fixed=T)
binarysites_full_max2[,2:57] <- sapply(binarysites_full_max2[2:57],as.numeric)
# Save this file; will be input for ASVs by site and marker
write.csv(binarysites_full_max2,"teleo_asvs_persample.csv",
          row.names=FALSE)

# Aggregate samples by site
aggsite_full_max <- binarysites_full_max2 %>% group_by(site) %>% 
  summarise_each(funs(sum(., na.rm=TRUE))) #this sums pres/abs counts for each taxon

# a) if wanting to look at detections across field replicates:
# Manually select which columns are fish
onlyfish_sep_max <- aggsite_full_max[,c(1,2,3,4,6,8,10,12,13,14,18,21,23,24,
                                26:29,33,35,36,38,
                                41:43,47,51:53,57)]
# Gather 
library(tidyr)
onlyfish_sep_max2 <- gather(onlyfish_sep_max, taxon, pres_abs, 
                        c("Salvelinus":"Salmo trutta"), 
                        factor_key=TRUE)

# Formatting column names to combine multiple ASVs per species:
# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
onlyfish_sep_max2$taxon <- sub("[0-9]$","",onlyfish_sep_max2$taxon)
# Get rid of spaces at the end of taxa
onlyfish_sep_max2$taxon <- sub("[ ]$","",onlyfish_sep_max2$taxon)
# Merge by taxa and add pres_abs
onlyfish_sep_max2$site <- as.factor(onlyfish_sep_max2$site)
onlyfish_sep_max2$taxon <- as.factor(onlyfish_sep_max2$taxon)
onlyfish_sep_max2$pres_abs <- as.numeric(onlyfish_sep_max2$pres_abs)
onlyfish_sep_max2<- onlyfish_sep_max2 %>% group_by(site,taxon) %>% 
  summarise_each(funs(sum))
# Save as input for work classifying certainty levels based on number of field detections
write.csv(onlyfish_sep_max2,"onlyfish_sep_max2.csv",row.names=FALSE)

# b) If wanting to look at just pres/abs across all field reps per site:
# Convert back to just pres/abs (i.e. pres if >0, abs if 0)
binaryaggsite_full_max <-as.data.frame(cbind(aggsite_full_max$site,(+(aggsite_full_max[-1] > 0))))
colnames(binaryaggsite_full_max)[1] <- "site"
# Manually select which columns are fish
onlyfish_full_max <- binaryaggsite_full_max[,c(1,2,3,4,6,8,10,12,13,14,18,21,23,24,
                                       26:29,33,35,36,38,
                                       41:43,47,51:53,57)]

# Formatting column names to combine multiple ASVs per species:
library(tidyr)
# Gather 
onlyfish_full_max2 <- gather(onlyfish_full_max, taxon, pres_abs, 
                         c("Salvelinus":"Salmo trutta"), 
                         factor_key=TRUE)
# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
onlyfish_full_max2$taxon <- sub("[0-9]$","",onlyfish_full_max2$taxon)
# Get rid of spaces at the end of taxa
onlyfish_full_max2$taxon <- sub("[ ]$","",onlyfish_full_max2$taxon)
# Merge by taxa and add pres_abs
onlyfish_full_max2$site <- as.factor(onlyfish_full_max2$site)
onlyfish_full_max2$taxon <- as.factor(onlyfish_full_max2$taxon)
onlyfish_full_max2$pres_abs <- as.numeric(onlyfish_full_max2$pres_abs)
uniqASVfish_max<- onlyfish_full_max2 %>% group_by(site,taxon) %>% 
  summarise_each(funs(sum))
# Spread back out again
uniqASVfish_max <- spread(uniqASVfish_max, key = taxon, value = pres_abs)
# This is the point where the 12Steleo data can be combined with other marker data(???)
write.csv(uniqASVfish_max,"uniqASVfish_maxlab_12steleo.csv",row.names=FALSE)


#-----Read counts/taxon/field sample AND total raw read depth/field sample---------
# drop stub column
lab_blank_max_sub <- lab_blank_max_sub[,-1]
# filter so that singletons are discarded (i.e. presence= read count >1)
lab_blank_max_sub[,2:57][lab_blank_max_sub[,2:57] <2] <- 0
# Convert NAs to zeros as well
lab_blank_max_sub[,2:57][is.na(lab_blank_max_sub[,2:57])] <- 0

# Get total raw read depth per sample:
asv_12s_full <- read.csv(here("Compute Canada Transfers","12Steleo","asv_counts_12steleo_full.csv"))
colnames(asv_12s_full)[1] <- "site"
# clean up site names
asv_12s_full$site <- gsub("-12Steleo_merged.assembled", "", asv_12s_full$site)
asv_12s_full$site <- gsub("12S-", "", asv_12s_full$site)
asv_12s_full$site <- gsub("-", "", asv_12s_full$site)

asv_12s_full$tot_depth <- rowSums(asv_12s_full[,2:7869])
#Cut to just total depth
sample_depth <- asv_12s_full[,c(1,7870)]

# Merge sample_depth df with lab_blank_sub_max df
teleo_depths <- merge(lab_blank_max_sub,sample_depth,by=c("site"),
                      all.x=TRUE,all.y=FALSE)

# Keep only fish
teleo_fish_depths <- teleo_depths[,c(1,2,3,4,6,8,10,12,13,14,18,21,23,24,
                                     26:29,33,35,36,38,
                                     41:43,47,51:53,57,58)]
colnames(teleo_fish_depths) <- gsub(x = colnames(teleo_fish_depths), #remove periods from column names
                                    pattern = ".", replacement = " ",fixed=T)


# Formatting column names to combine multiple ASVs per species:
# Gather 
library(tidyr)
teleo_fish_depths2 <- gather(teleo_fish_depths, taxon, depth, 
                             c("Salvelinus":"Salmo trutta"), 
                             factor_key=TRUE)

# Get rid of numbers at the end of taxa (i.e. "Aguilla 1")
teleo_fish_depths2$taxon <- sub("[0-9]$","",teleo_fish_depths2$taxon)
# Get rid of spaces at the end of taxa
teleo_fish_depths2$taxon <- sub("[ ]$","",teleo_fish_depths2$taxon)
# Merge by taxa and add pres_abs
teleo_fish_depths2$site <- as.factor(teleo_fish_depths2$site)
teleo_fish_depths2$taxon <- as.factor(teleo_fish_depths2$taxon)
teleo_fish_depths2$depth <- as.numeric(teleo_fish_depths2$depth)
teleo_fish_depths3 <- teleo_fish_depths2 %>% group_by(site,taxon,tot_depth) %>% 
  summarise_each(funs(sum))
teleo_fish_depths3$marker <- "12Steleo"

write.csv(teleo_fish_depths3,"12Steleo_depths_perfieldsample.csv",row.names=FALSE)

