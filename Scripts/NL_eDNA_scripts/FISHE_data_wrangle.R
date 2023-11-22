# This script is for wrangling the large FISHE asv counts output file to something more manageable 
# Use BLAST output for vertebrates and select only those ASVs that were >=99% SSS score
# Much more manageable to work with/load onto personal computer

library(readr)
library(tidyr)
blastfish_full <- read_delim("/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/blast_FISHE_range.out", delim="\t", col_names=F))
blastfish_full <- separate(blastfish_full, X1, into=c("readID","depth_allsamp"), sep=";")

# Create column for "sequence similarity score (sss) based on percent identity (pident) and percent query coverage (qcovs)
colnames(blastfish_full) <- c("readID", "depth_allsamp", "sacc", "ssciname", "sallseqid", "evalue","bitscore", "pident","qcovs")
blastfish_full$sss <- (blastfish_full$pident)*(blastfish_full$qcovs)/100

# Create new data frame only selecting rows with >=99% sss
library(dplyr)
speciesfish_full <- filter(blastfish_full,sss>=99)

# Create data frame for each ASV with sss>=99 and number of different taxonomic hits ("count")
asv_list <- data.frame(speciesfish_full %>% group_by(readID) %>%
                                         summarise(count = n_distinct(ssciname)))

# Load in ASV table:
counts <- read_tsv("/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/asv_counts_FISHE_range.tsv")

colnames(counts)[1] <- "readID"
merged <- merge(asv_list,counts, all.x=TRUE, all.y=FALSE)

# Save reduced ASV table file:
write.csv(merged, "reduced_asv_counts.csv", row.names=FALSE)

# Additionally, wrangle 12Steleo and MIFISHU ASV tables as well (convert to .csv and add column names etc.)
blast12Steleo_full <- read_tsv("/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/blast_12Steleo_full.out")
blast12Steleo_full <- separate(blast12Steleo_full, X1, into=c("readID","depth_allsamp"), sep=";")
colnames(blast12Steleo_full) <- c("readID", "depth_allsamp", "sacc", "ssciname", "sallseqid", "evalue","bitscore", "pident","qcovs")
write.csv(blast12Steleo_full, "blast_12Steleo_full.csv", row.names=F)

blastMIFISH_full <- read_tsv("/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/blast_MIFISHU_full.out")
blastMIFISH_full <- separate(blastMIFISH_full, X1, into=c("readID","depth_allsamp"), sep=";")
colnames(blastMIFISH_full) <- c("readID", "depth_allsamp", "sacc", "ssciname", "sallseqid", "evalue","bitscore", "pident","qcovs")
write.csv(blastMIFISH_full, "blast_MIFISHU_full.csv", row.names=F)
