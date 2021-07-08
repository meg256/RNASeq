#!/bin/bash

#Downloads and trims SRA reads. Runs FastQC before and after trimming.
#requires samples.txt as input - a line-separated list of sample names that SRA tools will fetch

##############################
#GET READS AND CONVERT
##############################
mkdir sra
cd sra
#make or move samples.txt into sra directory
while read line; do prefetch $line; done<samples.txt

#convert to fastq.gz
mkdir fastq
for file in *.sra; do parallel-fastq-dump --gzip --threads 8 --outdir ./fastq/ -s $file; done

#delete sra files if desired
#rm *.sra

#PRE-TRIM QC
#start in fastq folder
cd fastq
for file in *.fastq.gz;do fastqc $file;done
mkdir preQC
mv *.zip ./preQC
mv *.html ./preQC

#TRIM
#average four bases at a time over window. Trim at spots where average<20
#drop any reads that are less than 75bp
for file in *fastq.gz;do trimmomatic SE $file ${base}_TRIMMED.fastq.gz -threads 8 SLIDINGWINDOW:4:20 MINLEN:75
mkdir trimmed
mv *_TRIMMED.fastq.gz ./trimmed

#POST-TRIM QC
#start in fastq folder
cd trimmed
for file in *.fastq.gz;do fastqc $file;done
mkdir postQC
mv *.zip ./postQC
mv *.html ./postQC

