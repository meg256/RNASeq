library(DESeq2)
library(ggplot2)
counts <- as.matrix(read.csv("gene_count.txt",sep="\t",row.names=1,header=FALSE))

#read in coding for samples. Note that 1=control HF diet, 0=tumor HF diet
coldata <- read.csv("samplecode.txt",sep="\t",row.names=1,header=TRUE)

colnames(counts) <- rownames(coldata)

dds <- DESeqDataSetFromMatrix(countData = counts, colData, coldata, design = ~ Code)

#######################################
# PCA for 100 most variable genes (remove ntop=100 for default 500 most variable)
########################################


vsd <- vst(dds,blind=FALSE)

pcaData <- plotPCA(vsd, intgroup=c("Code"),ntop=100, returnData=TRUE)

percentVar <- round(100*attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1,PC2,color=Code)) + geom_point(size=3) + xlim(-2.5,2.5) + ylim(-1,1) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + geom_text(aes(label=name),vjust=2)

ggsave("PCA_100.png")



##########################################
#DE ANALYSIS
###########################################
dds<-DESeq(dds)

dds$Code <- factor(dds$Code, levels=c("wt","mu"))
                   
write.table(counts(dds,normalized=TRUE),file='Part1_normCounts.txt',sep="\t",row.names=TRUE,col.names=TRUE)
 
