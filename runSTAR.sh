#!/bin/bash


#preload ref genome into database
#STAR --genomeLoad LoadAndExit --genomeDir ./sra/ref
#would then need to add the following to the quant command:
#--genomeLoad LoadAndKeep --limitBAMsortRAM 4000000000 \ #limit of 4GB for BAM sorting step

#index ref genome
STAR --runMode genomeGenerate \
--runThreadN2 \
--genomeDir ./sra/ref \
--genomeFastaFiles ref.fa \ 
--sjdbGTFfile ref.gtf \
--sjdbOverhang 49

#NOTE: this loads UNTRIMMED reads. To map trimmed reads you must change filepath

#STAR map for SRS1794108
STAR --quantMode
--genomeDir ./sra/ref\
--runThreadN 2 \
--readFilesIn ./sra/fastq/SRS1794108.fastq.gz \ 
--readFilesCommand zcat \
--outFileNamePrefix control_SRS1794108 \ 
--outFilterMultimapNmax 1 \
--outReadsUnmapped unmapped_control \
--outSAMtype BAM SortedByCoordinate

#STAR map for SRS1794110
STAR --quantMode
--genomeDir ./sra/ref\
--runThreadN 2 \
--readFilesIn ./sra/fastq/SRS1794110.fastq.gz \ 
--readFilesCommand zcat \
--outFileNamePrefix control_SRS1794110 \ 
--outFilterMultimapNmax 1 \
--outReadsUnmapped unmapped_control_SRS1794110 \
--outSAMtype BAM SortedByCoordinate

#STAR map for SRS1794106
STAR --quantMode
--genomeDir ./sra/ref\
--runThreadN 2 \
--readFilesIn ./sra/fastq/SRS1794106.fastq.gz \ 
--readFilesCommand zcat \
--outFileNamePrefix control_SRS1794106 \ 
--outFilterMultimapNmax 1 \
--outReadsUnmapped unmapped_control_SRS1794106 \
--outSAMtype BAM SortedByCoordinate

#STAR map for SRS1794105
STAR --quantMode
--genomeDir ./sra/ref\
--runThreadN 2 \
--readFilesIn ./sra/fastq/SRS1794105.fastq.gz \ 
--readFilesCommand zcat \
--outFileNamePrefix tumor_SRS1794105 \ 
--outFilterMultimapNmax 1 \
--outReadsUnmapped unmapped_tumor_SRS1794105 \
--outSAMtype BAM SortedByCoordinate

#STAR map for SRS1794101
STAR --quantMode
--genomeDir ./sra/ref\
--runThreadN 2 \
--readFilesIn ./sra/fastq/SRS1794101.fastq.gz \ 
--readFilesCommand zcat \
--outFileNamePrefix tumor_SRS1794101 \ 
--outFilterMultimapNmax 1 \
--outReadsUnmapped unmapped_tumor_SRS1794101 \
--outSAMtype BAM SortedByCoordinate

#STAR map for SRS1794111
STAR --quantMode
--genomeDir ./sra/ref\
--runThreadN 2 \
--readFilesIn ./sra/fastq/SRS1794111.fastq.gz \ 
--readFilesCommand zcat \
--outFileNamePrefix tumor_SRS1794111_ \ 
--outFilterMultimapNmax 1 \
--outReadsUnmapped unmapped_tumor_SRS1794111 \
--outSAMtype BAM SortedByCoordinate

#OUTPUT for each - Log.final.out and ReadsPerGene.out.tab
#col 1 = gene ID, col 2 = counts for unstranded RNAseq
#col3 = counts for + strand RNA (for 3'RNA seq data)
#col4 = counts for - strand RNA (for stranded RNA seq data)