#!/bin/bash

cd ./sra/fastq/trimmed

#index BAM file with samtools
for f in *_Aligned.sortedByCoord.out.bam; do samtools index $f;done
#creates .bai file for each .bam


#creates tab delimited parallel merged mega-matrix for raw counts
paste *_ReadsPerGene.out.tab | \
cut -f1,2,6,10,14,18,22
tail -n +5 > gene_count.txt