#!/bin/bash
set -euo pipefail

# ========================================================================
# Advanced Bioinformatics Course Assignment 2026
# Student ID: SGUL: 2506253 KCL: 25163711
# =========================================================================

# MODIFIED VARIANT CALLER NGS PIPELINE 

# =========================================================================
# TOOLS AND DEPENDENCIES
# If not previously installed run once before beginning the pipeline
# ==========================================================================

# --- 1. INSTALLING ANACONDA ---

#wget https://repo.anaconda.com/archive/Anaconda3-2022.10-Linux-x86_64.sh
#chmod +x ./Anaconda3-2022.10-Linux-x86_64.sh
#bash ./Anaconda3-2022.10-Linux-x86_64.sh
#source ~/.bashrc

# --- 2. ADDING CONDA CHANNELS ---

#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge

# --- 3. NECESSARY BIOINFORMATIC TOOLS ---

#conda install samtools
#conda install bwa
#conda install freebayes
#conda install picard
#conda install bedtools
#conda install trimmomatic
#conda install fastqc
#sudo apt install libvcflib-tools

# --- 4. BCFTOOLS --- 
# Replaces freebayes in this modified pipeline and requires the following installations.
#conda install bcftools
#conda install -c conda-forge libopenblas

# --- 5. ANNOVAR ---

# Requires manual download from: https://www.openbioinformatics.org/annovar/annovar_download_form.php
# Downloaded link can be transfered to home directory via FileZilla. 
# Once downloaded, file requires extraction: 
#tar -zxvf ~/annovar.latest.tar.gz

# After extraction ANNOVAR annotation databases can be set up
#cd ~/annovar
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
#./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

# --- 6. SnpEff & SnpSift ---

# SnpEff Version 4.3 is required (as version 5.1 is produced errors with Java 8)
#conda install -c bioconda snpeff=4.3

#SnpEff also requires a hg19 annotation database to be downloaded
#java -jar ~/anaconda3/share/snpeff-4.3.1t-0/snpEff.jar download hg19

# Alongside SnpEff version 4.3, Snpsift version 4.3 is also required. 
# conda install -c bioconda snpsift=4.3

# =================================================================================
# MAKING PIPELINE DIRECTORY STRUCTURE
# Enables data organisation and efficient storage for data output
# =================================================================================

mkdir -p ~/ngs_pipeline
mkdir -p ~/ngs_pipeline/dnaseq
mkdir -p ~/ngs_pipeline/dnaseq/data
mkdir -p ~/ngs_pipeline/dnaseq/meta
mkdir -p ~/ngs_pipeline/dnaseq/results
mkdir -p ~/ngs_pipeline/dnaseq/logs
mkdir -p ~/ngs_pipeline/dnaseq/data/untrimmed_fastq
mkdir -p ~/ngs_pipeline/dnaseq/data/trimmed_fastq
mkdir -p ~/ngs_pipeline/dnaseq/data/reference
mkdir -p ~/ngs_pipeline/dnaseq/data/aligned_data
mkdir -p ~/ngs_pipeline/dnaseq/results/fastqc_untrimmed_reads
mkdir -p ~/ngs_pipeline/dnaseq/results/fastqc_trimmed_reads

# ==================================================================================
# NGS PIPELINE BEGINS 
# ==================================================================================

echo "This NGS pipeline runs FASTQC, read trimming, read alignment, variant calling, and variant annotation on raw sequencing data."

# --- 1. Download Data ---

# Download the initial input data files.
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed

# Rename the input FASTQ files from .qz to .gz extension
mv NGS0001.R1.fastq.qz NGS0001.R1.fastq.gz
mv NGS0001.R2.fastq.qz NGS0001.R2.fastq.gz

# Move input data files into correct directionarys
mv *.fastq.gz ~/ngs_pipeline/dnaseq/data/untrimmed_fastq
mv annotation.bed ~/ngs_pipeline/dnaseq/data

# Defining the raw FASTQ input file paths for quality control and trimming.
raw_R1=~/ngs_pipeline/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz
raw_R2=~/ngs_pipeline/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz

# Download the reference genome FASTA file and move it to the correct directory
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz
mv hg19.fa.gz ~/ngs_pipeline/dnaseq/data/reference

# Index of reference genome
# The index of the reference genome must be generated to enable alignment
# NOTE: bwa index takes  approx. 45 mins to run
bwa index ~/ngs_pipeline/dnaseq/data/reference/hg19.fa.gz

# --- 2. FASTQC on Untrimmed Data ---

# Run FASTQC to produce sequence quality reports on the untrimmed data
fastqc -t 4 $raw_R1 $raw_R2

# Move FASTQC output files into the untrimmed results directory
mv ~/ngs_pipeline/dnaseq/data/untrimmed_fastq/*fastqc* ~/ngs_pipeline/dnaseq/results/fastqc_untrimmed_reads/

# --- 3. Trimming Reads ---

# Use Trimmomatic tool to trim away adaptors and filter out reads.

trimmomatic PE \
  -threads 4 \
  -phred33  \
  $raw_R1 $raw_R2 \
  -baseout ~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.40-hdfd78af_0/share/trimmomatic-0.40-0/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

# Remove untrimmed fastq data to enable more storage space.
rm ~/ngs_pipeline/dnaseq/data/untrimmed_fastq/*.fastq.gz

# --- 4. FASTQC on Trimmed Data ---

# Run FASTQC on trimmed paired-end files
fastqc -t 4 \
~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P ~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P

# Move FastQC output into the trimmed results directory
mv ~/ngs_pipeline/dnaseq/data/trimmed_fastq/*fastqc* ~/ngs_pipeline/dnaseq/results/fastqc_trimmed_reads/

# --- 5. Read Alignment ---

# Note: The BWA index of hg19.fa.gz previously generated will be used in this step.

# Align trimmed paired end reads to hg19 reference genome, using BWA-MEM and read group information obtained from the FASTQ file header.
bwa mem -t 7 -v 1 \
  -R '@RG\tID:11V6WR1.111.D1375ACXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tPU:D1375ACXX.1' \
-I 250,50 \
~/ngs_pipeline/dnaseq/data/reference/hg19.fa.gz \
~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_1P \
~/ngs_pipeline/dnaseq/data/trimmed_fastq/NGS0001_trimmed_R_2P \
> ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001.sam

# Remove trimmed fastq data as it is no longer needed and enables more storage
rm ~/ngs_pipeline/dnaseq/data/trimmed_fastq/*

# --- 6. SAM to BAM Conversion, Sorting and Indexing ---

# Convert SAM file to BAM format
samtools view -h -b ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001.sam > ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001.bam

# Remove SAM file to enable more storage space
rm ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001.sam

# Sort the BAM file
samtools sort ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001.bam > ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted.bam

# Index the sorted bam file to generate a bai index file
samtools index ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted.bam

# --- 7. Duplicate Marking --- 

# Mark duplicate reads in sorted BAM
picard MarkDuplicates \
I=~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted.bam \
O=~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked.bam \
M=~/ngs_pipeline/dnaseq/data/aligned_data/marked_dup_metrics.txt

# Index the duplicate marked BAM file
samtools index ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked.bam

# Remove unneeded BAM and index files to enable more storage space.
rm ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted.bam
rm ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted.bam.bai

# --- 8. Post-Alignment Read Filtering ---

# Filter duplicate marked BAM file based on mapping quality and bitwise flags using samtools
samtools view -F 1796  -q 20 -o ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam \
~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked.bam

# Generate Index
samtools index ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam

# Remove unneeded BAM and index files to enable more storage space.
rm ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked.bam
rm ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked.bam.bai

# ---9. Standard Alignment Statistics --- 

# Summarise aligned read counts
samtools flagstat ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam \
> ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_flagstats.txt

# Generate mapping statistics of each chromosome
samtools idxstats ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam \
> ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_idxstats.txt

# Calculate insert size distribution
picard CollectInsertSizeMetrics \
  I=~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam \
  O=~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_insert_size_metrics.txt \
  H=~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_insert_size_histogram.pdf

# Calculate depth coverage 
samtools depth -b ~/ngs_pipeline/dnaseq/data/annotation.bed \
~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam \
> ~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_depth.txt

# --- 10. Variant Calling ---

# NOTE: Bcftools and libopenblas will have needed to be installed for this section to run.

# Decompress hg19 reference genome
zcat ~/ngs_pipeline/dnaseq/data/reference/hg19.fa.gz > ~/ngs_pipeline/dnaseq/data/reference/hg19.fa

# Index uncompressed reference genome using samtools faidx
samtools faidx ~/ngs_pipeline/dnaseq/data/reference/hg19.fa

# Generate genotype likelihoods with bcftools mpileup and call variants with bcftools call
bcftools mpileup -Ou -f ~/ngs_pipeline/dnaseq/data/reference/hg19.fa \
~/ngs_pipeline/dnaseq/data/aligned_data/NGS0001_sorted_marked_filtered.bam \
  | bcftools call -mv -Ov -o ~/ngs_pipeline/dnaseq/results/NGS0001_bcftools.vcf

# Compress and index VCF file
bgzip ~/ngs_pipeline/dnaseq/results/NGS0001_bcftools.vcf
tabix -p vcf ~/ngs_pipeline/dnaseq/results/NGS0001_bcftools.vcf.gz

# --- 11. Variant Filtering ---

# Filter low quality variants based on quality score, read depth and mapping quality.
vcffilter -f "QUAL > 20 & DP > 10 & MQ > 30" ~/ngs_pipeline/dnaseq/results/NGS0001_bcftools.vcf.gz > ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bcftools.vcf

# 11.1 Intersect VCF with BED file

# Using bedtools intersect to only retain variants overlapping target regions within annotation.bed.
bedtools intersect -header -wa -a ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bcftools.vcf \
-b ~/ngs_pipeline/dnaseq/data/annotation.bed \
> ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.vcf

# Compress and index filtered VCF file
bgzip ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.vcf
tabix -p vcf ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.vcf.gz

# Remove the uncompressed reference genome file to enable more storage space.
rm ~/ngs_pipeline/dnaseq/data/reference/hg19.fa

# --- 12. Variant Annotation ---

# 12.1 ANNOVAR
# Convert the filtered VCF file to .avinput format, compatible with ANNOVAR
~/annovar/convert2annovar.pl -format vcf4 \
~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.vcf.gz \
> ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.avinput

# Run ANNOVAR table function to annotate variants against multiple databases and output a csv file.
~/annovar/table_annovar.pl \
~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.avinput ~/annovar/humandb/ -buildver hg19 \
-out ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools -remove \
-protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout

# 12.2 SnpEff

# NOTE: Should have SnpEFF version 4.3 and snpEFF hg19 annotation database downloaded
# Annotate variants with predicted functional effects using snpEff hg19 database
java -Xmx8g -jar ~/anaconda3/share/snpeff-4.3.1t-0/snpEff.jar hg19 \
-stats ~/ngs_pipeline/dnaseq/results/NGS0001_snpeff_summary.html \
~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools.vcf.gz \
> ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools_annotated_snpeff.vcf

# --- 13. Variant Priotisation ---

# 13.1 SnpSift

# NOTE: SnpSift version 4.3 is required for this step to run.
# Performing basic variant priotisation on annotated SnpEff output using SnpSift. 
# Filter for exonic variants not in dbSNP 

java -Xmx8g -jar ~/anaconda3/share/snpsift-4.3.1t-3/SnpSift.jar filter \
"((ANN[*].IMPACT has 'HIGH') | (ANN[*].IMPACT has 'MODERATE')) & !(ID =~ 'rs')" \
~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools_annotated_snpeff.vcf \
> ~/ngs_pipeline/dnaseq/results/NGS0001_filtered_bedfile_bcftools_snpsift_prioritised.vcf

echo "NGS PIPELINE ENDED."
