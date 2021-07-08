# RNASeq
A workflow for processing SE RNA-seq data, from raw reads to DE visualization.
With credit/gratitude to the Griffith lab, Bioconductor RNA-seq workflow, the Cornell BioHPC Workshop for RNA-seq data analysis, and the wonderfully informative crazyhottommy walkthroughs (https://github.com/crazyhottommy/RNA-seq-analysis).

<b>As of July 2021:</b>
get_trim_QC.sh: downloads and trims SRA reads, runs FastQC before and after trimming. 
runSTAR.sh: preloads reference genome, then runs STAR for all samples in parallel. Outputs raw counts matrix as gene_counts.txt.  
cat.sh: Indexes BAM files as BAI. Combines raw count matrices into one big matrix
DE_PCA.R: takes in raw count matrix. Runs DE analysis. Outputs PCA and DE table.

<b>Practice dataset: </b>
SRA PRJNA353374
Penrose HM, Heller S, Cable C, Nakhoul H, Baddoo M, Flemington E, Crawford SE, Savkovic SD. High-fat diet induced leptin and Wnt expression: RNA-sequencing and pathway analysis of mouse colonic tissue and tumors. Carcinogenesis. 2017 Mar 1;38(3):302-311. doi: 10.1093/carcin/bgx001. PMID: 28426873; PMCID: PMC5862315. 
High-fat control (3 samples) and high-fat tumor (3 samples)

<b>Reference genome:</b>
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/
```
mkdir sra/ref
cd sra/ref
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.9_GRCm39/GCA_000001635.9_GRCm39_genomic.gtf.gz
```
Rename to ref.fa and ref.gtf; store in workdir/sra/ref
```
        workdir
           |
          sra
        /     \
      fastq    ref
      /    \
   preQC   trimmed
              |
            postQC
                 
```

<b>WORKFLOW</b>
<i>PROCESSING</i>
SRA-tools: fetch reads from NCBI SRA. https://github.com/ncbi/sra-tools
FastQC: QC for raw reads. https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
Trimmomatic: trim/filter raw reads. http://www.usadellab.org/cms/?page=trimmomatic
STAR: base-to-base alignment/mapping of reads to reference. https://github.com/alexdobin/STAR

<i>DIFFERENTIAL EXPRESSION/SPLICING</i>
DESeq2: differential expression analysis from raw read counts. https://bioconductor.org/packages/release/bioc/html/DESeq2.html
rMATS: differential alternative splicing. http://rnaseq-mats.sourceforge.net/rmatsdockerbeta/ <b>(in progress) </b>

<i>CLUSTERING AND FUNCTIONAL ANALYSIS</i>
topGO: https://bioconductor.org/packages/release/bioc/html/topGO.html <b>(in progress)</b>
clusterProfiler: https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html <b>(in progress)</b>

<i>VISUALIZATION</i>
Cytoscape: network visualization software. https://cytoscape.org <b>(in progress)</b>
bigPint: nifty visualization package for DE. https://lindsayrutter.github.io/bigPint <b>(in progress)</b>
  Rutter, L., Moran Lauter, A.N., Graham, M.A. et al. Visualization methods for differential expression analysis. BMC Bioinformatics 20, 458 (2019). https://doi.org/10.1186/s12859-019-2968-1
  Fig 2 (cluster comparison with parallel coordinate plots)
  Fig 9 (scatterplot + cluster comparison)
  Fig 10 (Litre plot)


