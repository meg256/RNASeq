install.packages('pathfindR')
BiocManager::install("KEGGREST")
BiocManager::install("biomaRt")


library(pathfindR)
library(KEGGREST)
library(biomaRt)
#KEGGREST for kegg sets
#pathFinder for mapping DE genes to kegg pathways

#non-human pathfindR vignette: https://cran.r-project.org/web/packages/pathfindR/vignettes/non_hs_analysis.html

#alt workflow: https://ucdavis-bioinformatics-training.github.io/2018-June-RNA-Seq-Workshop/friday/enrichment.html

depath='topTable.txt'


########################################################
#OBTAINING PIN AND GENE SETS
########################################################

#https://cran.r-project.org/web/packages/pathfindR/vignettes/obtain_data.html
gsets_list <- get_gene_sets_list(source = "KEGG",org_code = "mmu")

mmu_kegg_genes <- gsets_list$gene_sets
mmu_kegg_descriptions <- gsets_list$descriptions

## Save both as RDS files for later use
saveRDS(mmu_kegg_genes, "mmu_kegg_genes.RDS")
saveRDS(mmu_kegg_descriptions, "mmu_kegg_descriptions.RDS")

##load if needed
mmu_kegg_genes <- readRDS("mmu_kegg_genes.RDS")
mmu_kegg_descriptions <- readRDS("mmu_kegg_descriptions.RDS")

#######################################################
#GET PIN AND CONVERT 
#######################################################
url <- "https://stringdb-static.org/download/protein.links.v11.0/10090.protein.links.v11.0.txt.gz"
path2file <- file.path(tempdir(check = TRUE), "STRING.txt.gz")
download.file(url, path2file)
## read STRING pin file
mmu_string_df <- read.table(path2file, header = TRUE)

## filter using combined_score cut-off value of 800
mmu_string_df <- mmu_string_df[mmu_string_df$combined_score >= 800, ]

## fix ids
mmu_string_pin <- data.frame(Interactor_A = sub("^10090\\.", "", mmu_string_df$protein1),
                             Interactor_B = sub("^10090\\.", "", mmu_string_df$protein2))
head(mmu_string_pin, 2)
mmu_ensembl <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")

converted <- getBM(attributes = c("ensembl_peptide_id", "mgi_symbol"),
                   filters = "ensembl_peptide_id",
                   values = unique(unlist(mmu_string_pin)),
                   mart = mmu_ensembl)
mmu_string_pin$Interactor_A <- converted$mgi_symbol[match(mmu_string_pin$Interactor_A, converted$ensembl_peptide_id)]
mmu_string_pin$Interactor_B <- converted$mgi_symbol[match(mmu_string_pin$Interactor_B, converted$ensembl_peptide_id)]
mmu_string_pin <- mmu_string_pin[!is.na(mmu_string_pin$Interactor_A) & !is.na(mmu_string_pin$Interactor_B), ]
mmu_string_pin <- mmu_string_pin[mmu_string_pin$Interactor_A != "" & mmu_string_pin$Interactor_B != "", ]

head(mmu_string_pin, 2)

# remove self interactions
self_intr_cond <- mmu_string_pin$Interactor_A == mmu_string_pin$Interactor_B
mmu_string_pin <- mmu_string_pin[!self_intr_cond, ]

# remove duplicated inteactions (including symmetric ones)
mmu_string_pin <- unique(t(apply(mmu_string_pin, 1, sort))) # this will return a matrix object

mmu_string_pin <- data.frame(A = mmu_string_pin[, 1],
                             pp = "pp",
                             B = mmu_string_pin[, 2])

path2SIF <- file.path(tempdir(), "mmusculusPIN.sif")
write.table(mmu_string_pin,
            file = path2SIF,
            col.names = FALSE,
            row.names = FALSE,
            sep = "\t",
            quote = FALSE)
path2SIF <- normalizePath(path2SIF)

#############################################################
#READ IN DE TABLE AND SUBSET
#############################################################
#pathfindR filters bonferroni-adjusted enrichment p values 
df <- read.csv(depath,sep='\t',header=TRUE)
head(df,2)

#subset to columns w/equivalent of "Gene Symbols", "Change Values", and "p values"
deex <- subset(de,select=c(mgi_symbol,AveExpr,adj.P.Val))
head(deex)


#could limit even further to genes with >2fold change
de2fold <- subset(de,AveExpr>=2.0,select=c(mgi_symbol,AveExpr,adj.P.Val))
head(de2fold)


knitr::kable(head(deex))

###############################################################
# RUN PATHFINDR
###############################################################

#NOTE: Change pin_name_path to path2SIF if using a non-human, non-mouse file
# and change gene_sets to "Custom"
deout <- run_pathfindR(input = deex,
                                convert2alias = FALSE,
                                gene_sets = "mmu_KEGG",
                                #min_gset_size=5,max_gset_size=500,
                               # custom_genes = mmu_kegg_genes,
                               # custom_descriptions = mmu_kegg_descriptions,
                                pin_name_path = "mmu_STRING")

#default output of this function returns only terms containing at least 10 and at most 200 genes

#default subnetwork ssearch method: GREEDY with search depth of 1 and max depth of 1
# see https://github.com/egeulgen/pathfindR/wiki/Active-subnetwork-oriented-Enrichment-Documentation#selecting-the-active-subnetwork-search-algorithm 

#et voila, aggregated enrichment results for multiple filtered subnetworks
knitr::kable(deout)

################################################################
#VISUALIZATION
###############################################################
term_gene_graph(result_df=deout, use_description=TRUE)

#heatmap = can also convert to bar or boxplot
UpSet_plot(result_df=deout,genes_df=deex)
