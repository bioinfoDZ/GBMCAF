install.packages("XML")
BiocManager::install("graph")
BiocManager::install("GenomeInfoDb")
BiocManager::install("genefilter")
BiocManager::install("GSVA")
devtools::install_github("tidyverse/stringr")
install.packages('/Users/phillipgalbo/Desktop/GSVA_1.38.2.tar.gz', repos = NULL, lib = "/Library/Frameworks/R.framework/Versions/4.0/Resources/library")

library(genefilter)
library(GSVA)
library(Biobase)
library(reshape2)
library(stringr)

tcga_gbm_count <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/TCGA_GBM_counts.txt', header = T, row.names = 1, sep = '\t')
tcga_gbm_count <- 2^tcga_gbm_count
tcga_gbm_count[1:5,1:5]

gsva_gene_set <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/GSVA_Analysis/GSVA_geneset_AllCells.csv', header = T)
list <- split(as.matrix(gsva_gene_set)[,1], gsva_gene_set[,2])
tcga_gbm_gsva <- gsva(as.matrix(tcga_gbm), list, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = T) #Use if matrix is RPKM
tcga_gbm_gsva <- gsva(as.matrix(tcga_gbm_count), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T) #Use if matrix is Count
write.table(t(tcga_gbm_gsva), "/Users/phillipgalbo/Desktop/TCGA_GBM_gsva.txt", row.names=T, col.names = T, sep='\t')

tcga_lgg_count <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/TCGA_LGG_counts.txt', header = T, row.names = 1, sep = '\t')
tcga_lgg_count <- 2^tcga_lgg_count
tcga_lgg_count[1:5,1:5]
tcga_lgg_gsva <- gsva(as.matrix(tcga_lgg_count), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T) #Use if matrix is Count
write.table(t(tcga_lgg_gsva), "/Users/phillipgalbo/Desktop/TCGA_LGG_gsva.txt", row.names=T, col.names = T, sep='\t')

cgga_325 <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/CGGA.325.gene.count.name.txt', header = T, row.names = 1, sep = '\t')
cgga_325_gsva <- gsva(as.matrix(cgga_325), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T) #Use if matrix is Count
write.table(t(cgga_325_gsva), "/Users/phillipgalbo/Desktop/cgga_325_gsva.txt", row.names=T, col.names = T, sep='\t')

cgga_693 <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/CGGA.693.gene.count.name.txt', header = T, row.names = 1, sep = '\t')
cgga_693_gsva <- gsva(as.matrix(cgga_693), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T) #Use if matrix is Count
write.table(t(cgga_693_gsva), "/Users/phillipgalbo/Desktop/cgga_693_gsva.txt", row.names=T, col.names = T, sep='\t')

#Ivy atlas GSVA and data curation
gbm_ivyatlas <- read.csv(file.choose(), header = T, row.names = 1)
rows_genes <- read.csv(file.choose(), header = T, row.names = 1)
table(row.names(gbm_ivyatlas)==row.names(rows_genes))
row.names(gbm_ivyatlas) <- rows_genes$gene_symbol
write.table(gbm_ivyatlas, "/Users/phillipgalbo/Desktop/gbm_ivyatlas.txt", row.names=T, col.names = T, sep='\t')

gbm_ivyatlas <- read.table(file.choose(), header = T, row.names = 1, sep = '\t')
gbm_CAF_gene_set <- read.csv(file.choose(), header = T)
list <- split(as.matrix(gbm_CAF_gene_set)[,1], gbm_CAF_gene_set[,2])

gbm_ivyatlas_gsva <- gsva(as.matrix(gbm_ivyatlas), list, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = T) #Use if matrix is RPKM

gbm_ivyatlas_annotation <- read.csv(file.choose(),header = T)
gbm_ivyatlas_annotation$rna_well_id <- paste("X", gbm_ivyatlas_annotation$rna_well_id, sep= "")

gbm_ivyatlas_annotation_wanted <- gbm_ivyatlas_annotation[match(colnames(gbm_ivyatlas), gbm_ivyatlas_annotation$rna_well_id),]
dim(gbm_ivyatlas_annotation_wanted)

table(sort(colnames(gbm_ivyatlas_gsva)==gbm_ivyatlas_annotation_wanted$rna_well_id))

write.table(t(gbm_ivyatlas_gsva), "/Users/phillipgalbo/Desktop/gbm_ivyatlas_gsva.txt", row.names=T, col.names = T, sep='\t')
write.table(gbm_ivyatlas_annotation_wanted, "/Users/phillipgalbo/Desktop/gbm_ivyatlas_annotation_wanted.txt", row.names=T, col.names = T, sep='\t')



