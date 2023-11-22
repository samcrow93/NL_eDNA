#--------This script is for processing the NL eDNA metabarcoding data for 3 markers (12Steleo, MIFISHU, FISHE)------------------
# Originally run using Compute Canada

#load modules for cutadapt
module load StdEnv/2020 python/3.10.2
source /home/$USER/cutadapt_env/bin/activate
python -c "import cutadapt"

# 1.---- Primer removal----
#This will parallelize running through fastq files in a directory and demultiplex them by primer set (12Steleo, MIFISHU, or FISHE)
#Only use files that start with 12S (since FISHE and MIFISHU are triplicates within each site)
#fasta files with primer sequences for each marker, and params file, are in same directory
#creates new R1 and R2 fastq files that have been trimmed of their respective 5' primers (saves within subdirectory sep_markers/)
#discards pairs of reads if either of them do not contain primer sequences

# get parameters file
echo "projdir=/home/samcrow/scratch/eDNA/LAB_EDNA/" > cutadapt_params.tsv  # change this to whatever the path to your working directory is
source cutadapt_params.tsv

# change directory
cd $projdir

#export source variables
export projdir


#1. Primer removal
# First create file list of sample names (this will not include BLANK samples; they have a slightly different naming convention and will be processed in a different script)
mkdir sep_markers/
ls 12S*.fastq.gz | cut -d_ -f1 | sort | uniq > samples

# Then read in sample file, send to parallel (--joblog will create log file (cutadapt1.log) of which ones failed/succeeded)
cat samples | parallel --joblog cutadapt1.log ' cutadapt -g file:primers_fwd.Fasta -G file:primers_rev.Fasta -o sep_markers/{}_trim1-{name}_R1.fastq.gz -p sep_markers/{}_trim1-{name}_R2.fastq.gz --discard-untrimmed {}_R1.fastq.gz {}_R2.fastq.gz '

#1b) Blanks
#Blanks have slightly different naming convention than rest of files (have BLANK_number_), so cut first 14 character of sample name instead of after first _
ls 12S-BLANK*.fastq.gz | cut -b 1-14 | sort | uniq > blanks

cat blanks | parallel --jobs 32 ' cutadapt -g file:primers_fwd.Fasta -G file:primers_rev.Fasta -o sep_markers/{}_trim1-{name}_R1.fastq.gz -p sep_markers/{}_trim1-{name}_R2.fastq.gz --discard-untrimmed {}_R1.fastq.gz {}_R2.fastq.gz '


#then move these trimmed files into marker-specific subdirectories (MIFISHU, 12Steleo, FISHE)
cd sep_markers/
mkdir MIFISHU/
mkdir 12Steleo/
mkdir FISHE/
mv *trim1-MIFISHU*.fastq.gz MIFISHU/
mv *trim1-12Steleo*.fastq.gz 12Steleo/
mv *trim1-FISHE*.fastq.gz FISHE/

# 1.2 ----Remove reverse comp primers for 12Steleo data ----
# This code will parallelize running through samples for 12Steleo marker and trim primer sequence from the end of the read
# This is necessary because the 12Steleo fragment is short enough that sequencing runs through the read
# -m specifies minimum length of read (so -m 10 will discard anything less than length 10)

# get parameters file
echo "projdir=/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/" > cutadapt_params2.tsv
source cutadapt_params2.tsv

# change directory
cd $projdir

#export source variables
export projdir

mkdir trim2/
ls *.fastq.gz | cut -d_ -f1 | sort | uniq > samples

cat samples | parallel --jobs 32 ' cutadapt -a "CATGGTAAGTGTACCGGAAG; min_overlap=10" -A "AGAGTGACGGGCGGTGT; min_overlap=10" -m 10 -o trim2/{}_trim2-12Steleo_R1.fastq.gz -p trim2/{}_trim2-12Steleo_R2.fastq.gz {}_trim1-12Steleo_R1.fastq.gz {}_trim1-12Steleo_R2.fastq.gz '

# 1.2b) Blanks
ls 12S-BLANK*.fastq.gz | cut -b 1-14 | sort | uniq > blanks

cat blanks | parallel --jobs 16 ' cutadapt -a "CATGGTAAGTGTACCGGAAG; min_overlap=10" -A "AGAGTGACGGGCGGTGT; min_overlap=10" -m 10 -o trim2/{}_trim2-12Steleo_R1.fastq.gz -p trim2/{}_trim2-12Steleo_R2.fastq.gz {}_trim1-12Steleo_R1.fastq.gz {}_trim1-12Steleo_R2.fastq.gz '

#2. ------ Merging forward and reverse reads ------
# A) 12Steleo
# merging with PEAR; parallelize
# -q refers to minimum quality score; -v refers to minimum overlap (use 10 as minimum)
# places merged files in already created merged/ subdirectory (will create "assembled", "unassembled.forward", "unassembled.reverse", and "discarded" files

# get parameters file
echo "projdir=/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/" > 12steleo_merge_params.tsv
source 12steleo_merge_params.tsv

# change directory
cd $projdir

# export source variables
export projdir

# 1. Merge forward and reverse reads:
mkdir merged/
module load StdEnv/2020 pear/0.9.11
ls 12S-*.fastq.gz | cut -d_ -f1 | sort | uniq > tomerge

cat tomerge | parallel --jobs 32 --joblog 12s_merge.log ' pear -f {}_trim2-12Steleo_R1.fastq.gz -r {}_trim2-12Steleo_R2.fastq.gz -o /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/{}-12Steleo_merged -q 26 -v 10 '

# Merge reads for blanks
module load StdEnv/2020 pear/0.9.11
ls 12S-BLANK*fastq.gz | cut -b 1-14 | sort | uniq > blanktomerge

cat blanktomerge | parallel --jobs 32 ' pear -f {}_trim2-12Steleo_R1.fastq.gz -r {}_trim2-12Steleo_R2.fastq.gz -o /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/{}-12Steleo_merged -q 26 -v 10 '

# B) MIFISHU
# get parameters file
echo "projdir=/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/" > merge_params_mifishu.tsv
source merge_params_mifishu.tsv

# change directory
cd $projdir

# export source variables
export projdir

# 1. Merge forward and reverse reads:
mkdir merged/
module load StdEnv/2020 pear/0.9.11
ls 12S-*.fastq.gz | cut -d_ -f1 | sort | uniq > tomerge

cat tomerge | parallel --jobs 32 --joblog mifish_merge.log ' pear -f {}_trim1-MIFISHU_R1.fastq.gz -r {}_trim1-MIFISHU_R2.fastq.gz -o /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/{}-MIFISHU_merged -q 26 -v 10 '

# Merge reads for blanks
module load StdEnv/2020 pear/0.9.11
ls 12S-BLANK*fastq.gz | cut -b 1-14 | sort | uniq > blanktomerge

cat blanktomerge | parallel --jobs 32 ' pear -f {}_trim1-MIFISHU_R1.fastq.gz -r {}_trim1-MIFISHU_R2.fastq.gz -o /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/{}-MIFISHU_merged -q 26 -v 10 '

# C) FISHE
# get parameters file
echo "projdir=/home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/" > merge_params_fishe.tsv
source merge_params_fishe.tsv

# change directory
cd $projdir

# export source variables
export projdir

# 1. Merge forward and reverse reads:
mkdir merged/
module load StdEnv/2020 pear/0.9.11
ls 12S-*.fastq.gz | cut -d_ -f1 | sort | uniq > tomerge

cat tomerge | parallel --jobs 32 --joblog fishe_merge.log ' pear -f {}_trim1-FISHE_R1.fastq.gz -r {}_trim1-FISHE_R2.fastq.gz -o /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/{}-FISHE_merged -q 26 -v 10 '

# Merge reads for blanks
module load StdEnv/2020 pear/0.9.11
ls 12S-BLANK*fastq.gz | cut -b 1-14 | sort | uniq > blanktomerge

cat blanktomerge | parallel --jobs 32 ' pear -f {}_trim1-FISHE_R1.fastq.gz -r {}_trim1-FISHE_R2.fastq.gz -o /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/{}-FISHE_merged -q 26 -v 10 '

#3. --------------Data concatenation------------
# Combines data from each merged sample file into one fastq file
# A) 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
# Start by getting rid of unassembled and discarded files so they are not concatenated as well
rm *unassembled* *discarded*
# then concatenate assembled (merged) files:
for f in *;
do
sed -e "s/\(^@A00.*\) .*$/\1;sample=${f%.*};/" $f >> concat_12Steleo_full.fastq
done

# B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
rm *unassembled* *discarded*
for f in *;
do
sed -e "s/\(^@A00.*\) .*$/\1;sample=${f%.*};/" $f >> concat_MIFISHU_full.fastq
done

# C) FISHE
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
rm *unassembled* *discarded*
for f in *;
do
sed -e "s/\(^@A00.*\) .*$/\1;sample=${f%.*};/" $f >> concat_FISHE_full.fastq
done


# 4. -----------------Quality Filtering-------------------
# Use VSEARCH --fastx_filter
# parameters: -fastq_maxee refers to expected number of errors (set to 1, will exclude reads that have at least 1 expected error)
# A) 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
module load vsearch/2.21.1
vsearch --fastx_filter concat_12Steleo_full.fastq --fastq_maxee 1 --fastaout filt_12Steleo_full.fasta

# B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
vsearch --fastx_filter concat_MIFISHU_full.fastq --fastq_maxee 1 --fastaout filt_MIFISHU_full.fasta

# C) FISHE
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
vsearch --fastx_filter concat_FISHE_full.fastq --fastq_maxee 1 --fastaout filt_FISHE_full.fasta

# 5. --------------Dereplication-------------------------
# this step outputs the unique ASVs (--sizeout command records number of copies of each ASV in whole dataset)
#A) 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
vsearch --derep_fulllength filt_12Steleo_full.fasta --output derep_12Steleo_full.fasta --sizeout --relabel uniq

#B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
vsearch --derep_fulllength filt_MIFISHU_full.fasta --output derep_MIFISHU_full.fasta --sizeout --relabel uniq

#C) FISHE
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
vsearch --derep_fulllength filt_FISHE_full.fasta --output derep_FISHE_full.fasta --sizeout --relabel uniq

#6. ----------------Denoising----------------------------
# Uses read frequency and sequence composition to infer likely sequencing errors; compares ASVs against one another and applies frequency thresholds relative to counts
# this step also includes size filtering (size refering to number of reads per ASV; --minsize 8 is the default and will discard any ASVs that occur less than 8 times
#A) 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
vsearch --cluster_unoise derep_12Steleo_full.fasta --minsize 8 --unoise_alpha 2 --centroids denoise_12Steleo_full.fasta

#B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
vsearch --cluster_unoise derep_MIFISHU_full.fasta --minsize 8 --unoise_alpha 2 --centroids denoise_MIFISHU_full.fasta

#C) FISHE
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
vsearch --cluster_unoise derep_FISHE_full.fasta --minsize 8 --unoise_alpha 2 --centroids denoise_FISHE_full.fasta

# 7. --------------Length Filtering --------------------
#A) # 12Steleo: cut to  60-70 bp range
# have to unwrap data as output by vsearch first:
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' denoise_12Steleo_full.fasta > denoise_unwrap_12Steleo_full.fasta
# then filter so keep only reads within amplicon size range:
vsearch --fastx_filter denoise_unwrap_12Steleo_full.fasta --fastq_minlen 60 --fastq_maxlen 70 --fastaout size_filt_12Steleo_full.fasta

#B) MFISHU: cut to 163-185 bp range
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
perl -pe '$. > 1 and /^>/ ? print "\n" : chomp' denoise_MIFISHU_full.fasta > denoise_unwrap_MIFISHU_full.fasta
# then filter so keep only reads within amplicon size range:
vsearch --fastx_filter denoise_unwrap_MIFISHU_full.fasta --fastq_minlen 163 --fastq_maxlen 185 --fastaout size_filt_MIFISHU_full.fasta

#C) FISHE: cut to 223-232 bp range
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
vsearch --fastx_filter denoise_unwrap_FISHE_full.fasta --fastq_minlen 223 --fastq_maxlen 232 --fastaout size_filt_FISHE_range.fasta

# 8.---------------Chimaera filtering -----------------
# A). 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
vsearch --uchime3_denovo size_filt_12Steleo_full.fasta --nonchimeras nochim_12Steleo_full.fasta

#B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
vsearch --uchime3_denovo size_filt_MIFISHU_full.fasta --nonchimeras nochim_MIFISHU_full.fasta

#C) FISHE
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
vsearch --uchime3_denovo size_filt_FISHE_range.fasta --nonchimeras nochim_FISHE_range.fasta

#9. ----------------Generating ASV table-----------------
# this code searches all the reads in the concatenated file (not yet quality filtered or filtered in any way) and maps them to the ASVs as determined after the chimera filtering step
# first, convert concat.fastq file to fasta
#A. 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
module load fastx-toolkit/0.0.14
fastq_to_fasta -i concat_12Steleo_full.fastq -o concat_12Steleo_full.fasta
vsearch --search_exact concat_12Steleo_full.fasta -db nochim_12Steleo_full.fasta -otutabout asv_counts_12Steleo_full.tsv

#B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
fastq_to_fasta -i concat_MIFISHU_full.fastq -o concat_MIFISHU_full.fasta
vsearch --search_exact concat_MIFISHU_full.fasta -db nochim_MIFISHU_full.fasta -otutabout asv_counts_MIFISHU_full.tsv

#C) FISHE
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
fastq_to_fasta -i concat_FISHE_full.fastq -o concat_FISHE_full.fasta
vsearch --search_exact concat_FISHE_full.fasta -db nochim_FISHE_range.fasta -otutabout asv_counts_FISHE_range.tsv


#10. --------Taxonomic assignment of ASVs using BLAST---------------------
# Note that this uses BLAST nucleotide database available on Compute Canada (downloaded 23/03/2022)
#A). 12Steleo
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/12Steleo/trim2/merged/
module load StdEnv/2020 gcc/9.3.0 blast+/2.12.0
export BLASTDB=/cvmfs/bio-test.data.computecanada.ca/content/databases/Core/blast_dbs/2022_03_23

# -sorthits 3 specifies sorting hits by percent identity; -outfmt 6 specifies tabular output format
blastn -db nt -num_threads $SLURM_CPUS_PER_TASK -query nochim_12Steleo_full.fasta -taxidlist vertebrate.txids -sorthits 3 -perc_identity 95 -outfmt "6 qacc sacc ssciname sallseqid evalue bitscore pident qcovs" -out blast_12Steleo_full.out

#B) MIFISHU
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/MIFISHU/merged/
blastn -db nt -num_threads $SLURM_CPUS_PER_TASK -query nochim_MIFISHU_full.fasta -taxidlist vertebrate.txids -sorthits 3 -perc_identity 95 -outfmt "6 qacc sacc ssciname sallseqid evalue bitscore pident qcovs" -out blast_MIFISHU_full.out

#C) FISHE
# note that this marker picked up majority inverts; limit BLAST to just vertebrate taxonomy IDs
cd /home/samcrow/scratch/eDNA/LAB_EDNA/sep_markers/FISHE/merged/
# First, generate list of NCBI taxonomy IDs for just vertebrates:
sh -c "$(wget -q ftp://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
export PATH=${PATH}:${HOME}/edirect
esearch -db taxonomy -query "Vertebrata[ORGN]" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId > vertebrate.txids
# Then, run blast by limiting search to just these vertebrate taxids:
blastn -db nt -num_threads $SLURM_CPUS_PER_TASK -query nochim_FISHE_range.fasta -taxidlist vertebrate.txids -sorthits 3 -perc_identity 95 -outfmt "6 qacc sacc ssciname staxid sallseqid evalue bitscore pident qcovs" -out blast_FISHE_range.out
cp blast_FISHE_range.out blast_FISHE_range.csv

#  ----------------- Wrangling FISHE ASV table----------------------
# FISHE marker picked up majority invertebrates
# ASV table as output above is too unwieldy/large to open in Rstudio
# R script uses BLAST output for all taxa and selects only those ASVs taht were >=99% sequence similarity score (done in Rstudio for MIFISHU and 12Steleo, but need higher computing power for this marker)
# Resulting output is much more manageable to work with on personal computer
# Note that this script also formats blast.out files for MIFISHU and 12Steleo markers and converts them to .csv files
module load nixpkgs/16.09 gcc/7.3.0 StdEnv/2020 r/4.2.1
Rscript FISHE_data_wrangle.R


# At this point, ASV table files and blast output files can be transferred to personal computer
# remainder of analysis done in Rstudio







































































































