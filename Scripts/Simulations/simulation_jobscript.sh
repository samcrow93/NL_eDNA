#----------------------This script contains the code necessary for running metabarcoding sequencing simulations---------------
# Requires Grinder software (https://github.com/zyxue/biogrinder)

# 1)--------------Run Grinder for each marker--------------------------
# Simulate 4million reads using Illumina error model and at varying sequence abundances (see input *_abund_*.txt files)
# A. 12Steleo
# load required packages
module load gcc bioperl

# Run script from scratch/simulations/Grinder-0.5.4 directory, cd into simulations_abund directory (params file specifies this projdir)
# get parameters file
echo "projdir=/home/samcrow/scratch/simulations/Grinder-0.5.4/simulations_abund" > grinder_params_abund.tsv
source grinder_params_abund.tsv

# Change directory
cd $projdir

# export source variables
export projdir

# Run grinder to simulate 12Steleo amplicons:
# First, list file names:
ls lab*.fasta nfld*.fasta | cut -d. -f1 | sort | uniq > filenames_12Steleo_abund

cat filenames_12Steleo_abund | parallel --jobs 32 '../bin/grinder -rf {}.fasta -fr primers_12Steleo.fasta -length_bias 0 -tr 4000000 -af {}_abund_12Steleo.txt -id 150 -rd 150 -fq 1 -ql 36 30 -mo FR -dc '-' -un 1 -md poly4 3e-3 3.3e-8 -mr 98 2 -bn {}_out_12Steleo_ab '


# B. MIFISHU data
# Run grinder to simulate MIFISHU amplicons:
# First, list file names:
ls lab*.fasta nfld*.fasta | cut -d. -f1 | sort | uniq > filenames_mifish_abund

cat filenames_mifish_abund | parallel --jobs 32 '../bin/grinder -rf {}.fasta -fr primers_mifish.fasta -length_bias 0 -tr 4000000 -af {}_abund_mifish.txt -id 245 -rd 150 -fq 1 -ql 36 30 -mo FR -dc '-' -un 1 -md poly4 3e-3 3.3e-8 -mr 98 2 -bn {}_out_mifish_ab '


# C. FISHE data
# Run grinder to simulate FISHE amplicons:
# First, list file names:
ls lab*.fasta nfld*.fasta | cut -d. -f1 | sort | uniq > filenames_fishe_abund

cat filenames_fishe_abund | parallel --jobs 32 '../bin/grinder -rf {}.fasta -fr primers_fishe.fasta -length_bias 0 -tr 4000000 -af {}_abund_fishe.txt -rd 150 -id 281 -fq 1 -ql 36 30 -mo FR -dc '-' -un 1 -md poly4 3e-3 3.3e-8 -mr 98 2 -bn {}_out_fishe_ab '


# 1.5) ---------------Move output fastq files into separate directories---------------------
mkdir 12Steleo_ab
mkdir mifish_ab
mkdir fishe_ab
mv *_12Steleo_ab-reads.fastq 12Steleo_ab
mv *_mifish_ab-reads.fastq mifish_ab
mv *_fishe_ab-reads.fastq fishe_ab


# 2)--------- Separate interleaved R1 R2 fastq files output by Grinder -------------------
# A. 12Steleo
cd /home/samcrow/projects/def-ibradbur/samcrow/simulations/Grinder-0.5.4/simulations_abund/12Steleo_ab/
# Load required package
module load StdEnv/2020 seqtk/1.3

# get parameters file
echo "projdir=/home/samcrow/scratch/simulations/Grinder-0.5.4/simulations_abund_12Steleo_ab" > 12Steleo_separate_params_ab.tsv
source 12Steleo_separate_params_ab.tsv

# Change directory
cd $projdir

# Export source variables
export projdir

# Output of grinder is an interleaved fastq; need to separate R1 and R2 into separate files
# First, list file names:
ls *.fastq | cut -d_ -f1 | sort | uniq > files_to_sep_ab

cat files_to_sep_ab | parallel --jobs 32 ' seqtk seq -1 {}_out_12Steleo_ab-reads.fastq > {}_12Steleo_R1_ab.fastq '
cat files_to_sep_ab | parallel --jobs 32 ' seqtk seq -2 {}_out_12Steleo_ab-reads.fastq > {}_12Steleo_R2_ab.fastq '

#B. MIFISHU
# get parameters file
cd /home/samcrow/projects/def-ibradbur/samcrow/simulations/Grinder-0.5.4/simulations_abund/mifish_ab/
echo "projdir=/home/samcrow/scratch/simulations/Grinder-0.5.4/simulations_abund_mifish_ab" > mifish_separate_params_ab.tsv
source mifish_separate_params_ab.tsv

# Change directory
cd $projdir

# Export source variables
export projdir

# Output of grinder is an interleaved fastq; need to separate R1 and R2 into separate files
# First, list file names:
ls *.fastq | cut -d_ -f1 | sort | uniq > files_to_sep_mifish_ab

cat files_to_sep_mifish_ab | parallel --jobs 32 ' seqtk seq -1 {}_out_mifish_ab-reads.fastq > {}_mifish_R1_ab.fastq '
cat files_to_sep_mifish_ab | parallel --jobs 32 ' seqtk seq -2 {}_out_mifish_ab-reads.fastq > {}_mifish_R2_ab.fastq '

#C. FISHE
# get parameters file
cd /home/samcrow/projects/def-ibradbur/samcrow/simulations/Grinder-0.5.4/simulations_abund/fishe_ab/
echo "projdir=/home/samcrow/scratch/simulations/Grinder-0.5.4/simulations_abund/fishe_ab" > fishe_separate_params_ab.tsv
source fishe_separate_params_ab.tsv

# Change directory
cd $projdir

# Export source variables
export projdir

# Output of grinder is an interleaved fastq; need to separate R1 and R2 into separate files
# First, list file names:
ls *.fastq | cut -d_ -f1 | sort | uniq > files_to_sep_fishe_ab

cat files_to_sep_fishe_ab | parallel --jobs 32 ' seqtk seq -1 {}_out_fishe_ab-reads.fastq > {}_fishe_R1_ab.fastq '
cat files_to_sep_fishe_ab | parallel --jobs 32 ' seqtk seq -2 {}_out_fishe_ab-reads.fastq > {}_fishe_R2_ab.fastq '


# Following this, fastq files are ready to be processed with custom pipeline and DADA2 pipeline
# For custom pipeline, see metabarcode_master_script.sh in https://github.com/samcrow93/NL_eDNA/tree/main/Scripts/NL_eDNA_scripts
# For DADA2 pipeline, see https://benjjneb.github.io/dada2/bigdata_paired.html, with the following changes to parameters in filtering step:
# filterAndTrim(fwd=file.path(pathF, fastqFs), filt=file.path(filtpathF, fastqFs),
  #            rev=file.path(pathR, fastqRs), filt.rev=file.path(filtpathR, fastqRs),
   #           maxEE=2, minQ=2, truncQ=2, maxN=0, compress=FALSE, multithread=TRUE)













