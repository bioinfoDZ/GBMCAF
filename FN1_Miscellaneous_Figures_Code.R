#FN1 Missing Data Code
library(ggplot2)

#1: Get proteomic ranking
proteomic_ranking <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/FN1_Miscellaneous_Figures/Proteomic_Protein_Rank/proteomic_protein_rank.txt', header = T, sep = '\t') #read: proteomic_protein_rank.txt

ggplot(data=proteomic_ranking, aes(y=IMR90_CM_Average, x= IMR90_CM_Ranking)) + geom_point() + coord_cartesian(ylim = c(-5, 10)) + theme_classic() + geom_hline(yintercept=0) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))

ggplot(data=proteomic_ranking, aes(y=GBM_CAF_CM_Average, x= GBM_CAF_CM_Ranking)) + geom_point() + coord_cartesian(ylim = c(-5, 10)) + theme_classic() + geom_hline(yintercept=0) + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))

#2: Get FN1 gene expression correlation with CAF GSVA enrichment scores
tcga_gbm_fpkm <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/TCGA_GBM_RNAseq_no_duplicates_test.txt', header = T, row.names = 1, sep = '\t')
dim(tcga_gbm_fpkm)
tcga_gbm_gsva <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/GSVA_Analysis/TCGA/TCGA_GBM_COX_analysis.csv', header = T, row.names = 1)
tcga_gbm_fpkm_wanted <- tcga_gbm_fpkm[,match(row.names(tcga_gbm_gsva), colnames(tcga_gbm_fpkm))]
dim(tcga_gbm_fpkm_wanted)
secreted_factors <- c('FN1', 'CXCL8', 'CCL2', 'IL6')
tcga_gbm_fpkm_wanted_secreted_factors <- tcga_gbm_fpkm_wanted[match(secreted_factors, row.names(tcga_gbm_fpkm_wanted)),]
dim(tcga_gbm_fpkm_wanted_secreted_factors)
tcga_gbm_fpkm_wanted_secreted_factors[1:4,1:4]
tcga_gbm_fpkm_wanted_secreted_factors <- t(tcga_gbm_fpkm_wanted_secreted_factors)
table(sort(row.names(tcga_gbm_fpkm_wanted_secreted_factors)==row.names(tcga_gbm_gsva)))
tcga_gbm_fpkm_wanted_secreted_factors_correlation <- cbind(tcga_gbm_fpkm_wanted_secreted_factors, tcga_gbm_gsva$CAF)
name <- c('FN1', 'CXCL8', 'CCL2', 'IL6', 'CAF')
colnames(tcga_gbm_fpkm_wanted_secreted_factors_correlation) <- name
tcga_gbm_fpkm_wanted_secreted_factors_correlation <- data.frame(tcga_gbm_fpkm_wanted_secreted_factors_correlation)
ggplot(tcga_gbm_fpkm_wanted_secreted_factors_correlation, aes(x=FN1, y=CAF)) + 
  geom_point(shape=18, color="blue", size = 5)+
  geom_smooth(method=lm,  linetype="dashed", color="darkred", fill="blue") + labs(x = "FN1 Expression", y = "CAF Enrichment Score") + theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))
cor.test(tcga_gbm_fpkm_wanted_secreted_factors_correlation$FN1, tcga_gbm_fpkm_wanted_secreted_factors_correlation$CAF, method=c("pearson"))

cgga325 <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/GBM_CGGA325_Gene_Expression_Profile_FPKM.csv', header = T, row.names = 1)
dim(cgga325)
cgga_gbm_gsva <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/GSVA_Analysis/CGGA_325/CGGA_325_COX_analysis.csv', header = T, row.names = 1)
cgga_gbm_gsva_gbm <- subset(cgga_gbm_gsva, cgga_gbm_gsva$Grade == 4)
secreted_factors <- c('FN1', 'IL8', 'CCL2', 'IL6')
cgga325_secreted_factors <- cgga325[match(secreted_factors, row.names(cgga325)),]
cgga325_secreted_factors[1:4,1:4]
cgga325_secreted_factors <- t(cgga325_secreted_factors)
cgga325_secreted_factors <- cgga325_secreted_factors[match(row.names(cgga_gbm_gsva_gbm), row.names(cgga325_secreted_factors)),]
table(sort(row.names(cgga325_secreted_factors)==row.names(cgga_gbm_gsva_gbm)))
cgga325_secreted_factors_correlation <- cbind(cgga325_secreted_factors, cgga_gbm_gsva_gbm$CAF)
name <- c('FN1', 'CXCL8', 'CCL2', 'IL6', 'CAF')
colnames(cgga325_secreted_factors_correlation) <- name
cgga325_secreted_factors_correlation <- data.frame(cgga325_secreted_factors_correlation)
ggplot(cgga325_secreted_factors_correlation, aes(x=log2(FN1+1), y=CAF)) + 
  geom_point(shape=18, color="blue", size = 5)+
  geom_smooth(method=lm,  linetype="dashed", color="darkred", fill="blue") + labs(x = "FN1 Expression", y = "CAF Enrichment Score") + theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))
cor.test(log2(cgga325_secreted_factors_correlation$FN1+1), cgga325_secreted_factors_correlation$CAF, method=c("pearson"))

cgga693 <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/cgga_693_gbm_FPKM.txt', header = T, row.names = 1, sep = '\t')
dim(cgga693)
cgga693_gbm_gsva <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/GSVA_Analysis/CGGA_693/CGGA_693_COX_analysis.csv', header = T, row.names = 1)
cgga693_gbm_gsva <- subset(cgga693_gbm_gsva, cgga693_gbm_gsva$Grade == 'WHO IV')
secreted_factors <- c('FN1', 'IL8', 'CCL2', 'IL6')
cgga693_secreted_factors <- cgga693[match(secreted_factors, row.names(cgga693)),]
cgga693_secreted_factors[1:4,1:4]
cgga693_secreted_factors <- t(cgga693_secreted_factors)
cgga693_secreted_factors <- cgga693_secreted_factors[match(row.names(cgga693_gbm_gsva), row.names(cgga693_secreted_factors)),]
table(sort(row.names(cgga693_secreted_factors)==row.names(cgga693_gbm_gsva)))
cgga693_secreted_factors_correlation <- cbind(cgga693_secreted_factors, cgga693_gbm_gsva$CAF)
name <- c('FN1', 'CXCL8', 'CCL2', 'IL6', 'CAF')
colnames(cgga693_secreted_factors_correlation) <- name
cgga693_secreted_factors_correlation <- data.frame(cgga693_secreted_factors_correlation)
ggplot(cgga693_secreted_factors_correlation, aes(x=log2(FN1+1), y=CAF)) + 
  geom_point(shape=18, color="blue", size = 5)+
  geom_smooth(method=lm,  linetype="dashed", color="darkred", fill="blue") + labs(x = "FN1 Expression", y = "CAF Enrichment Score") + theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))
cor.test(log2(cgga693_secreted_factors_correlation$FN1+1), cgga693_secreted_factors_correlation$CAF, method=c("pearson"))

cptac_proteomic <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/CPTAC/cptac_proteome_normalized.txt', header = T, sep = '\t')
cptac_gsva <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/CPTAC/CPTAC_GSVA.txt', header = T, row.names = 1, sep = '\t')
cptac_proteomic_secreted_factors <- subset(cptac_proteomic, cptac_proteomic$symbol == 'FN1' | cptac_proteomic$symbol == 'CXCL8' | cptac_proteomic$symbol == 'CCL2' | cptac_proteomic$symbol == 'IL-6')
row.names(cptac_proteomic_secreted_factors) <- cptac_proteomic_secreted_factors$symbol
cptac_proteomic_secreted_factors <- cptac_proteomic_secreted_factors[,-1]
cptac_proteomic_secreted_factors <- t(cptac_proteomic_secreted_factors)
cptac_proteomic_secreted_factors <- cptac_proteomic_secreted_factors[match(row.names(cptac_gsva), row.names(cptac_proteomic_secreted_factors)),]
table(sort(row.names(cptac_proteomic_secreted_factors)==row.names(cptac_gsva)))
cptac_proteomic_secreted_factors_correlation <- cbind(cptac_proteomic_secreted_factors, cptac_gsva$CAF)
name <- c('CCL2', 'CXCL8', 'FN1', 'CAF')
colnames(cptac_proteomic_secreted_factors_correlation) <- name
cptac_proteomic_secreted_factors_correlation <- data.frame(cptac_proteomic_secreted_factors_correlation)
ggplot(cptac_proteomic_secreted_factors_correlation, aes(x=FN1, y=CAF)) + 
  geom_point(shape=18, color="blue", size=5)+
  geom_smooth(method=lm,  linetype="dashed", color="darkred", fill="blue") + labs(x = "FN1 Expression", y = "CAF Enrichment Score") + theme_classic() + theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"))

cor.test(cptac_proteomic_secreted_factors_correlation$FN1, cptac_proteomic_secreted_factors_correlation$CAF, method=c("pearson"))

#3: FN1 gene expression spatial transcriptomics plot
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)

gbm_spatial_data_0.6 <- readRDS(file.choose()) #read:gbm_spatial_data_0.6.rds
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "FN1")
VlnPlot(gbm_spatial_data_0.6, feature = 'FN1')

#4: Single Cell RNA-seq bubble plot (initial data matrix genreated on Deyou's computer)
gbm_seq <- read.table('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/normal_cells_with_cafs.txt', header = T, row.names = 1, sep = "\t")
sc <- readRDS('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/sc_rfcorrect.rds')
sc_umap <- readRDS( '/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/sc_umap.rds')

vascular_cells <- subset(sc_umap@umap, sc_umap@umap$V1 > -5 & sc_umap@umap$V1 < 5 & sc_umap@umap$V2 > -5 & sc_umap@umap$V2 < -2.5)
cluster_want_one <- subset(sc_umap@umap, sc_umap@umap$V1 > -3 & sc_umap@umap$V1 < 1 & sc_umap@umap$V2 > 0 & sc_umap@umap$V2 < 5)
cluster_want_two <- subset(sc_umap@umap, sc_umap@umap$V1 > -5 & sc_umap@umap$V1 < -2 & sc_umap@umap$V2 > 2 & sc_umap@umap$V2 < 5)

test <- sc_umap@cluster$kpart
test <- data.frame(test)
cluster_2_names <- subset(test, test ==2)
gbm_seq_cluster2 <- gbm_seq[,match(row.names(cluster_2_names), colnames(gbm_seq))]
dim(gbm_seq_cluster2)
gbm_seq_cluster2[1:5,1:5]

cluster5_names <- subset(test, test == 5)
gbm_seq_cluster5 <- gbm_seq [,match(row.names(cluster5_names), colnames(gbm_seq))]
dim(gbm_seq_cluster5)
gbm_seq_cluster5[1:5,1:5]

cluster6_names <- subset(test, test == 6)
gbm_seq_cluster6 <- gbm_seq [,match(row.names(cluster6_names), colnames(gbm_seq))]
dim(gbm_seq_cluster6)
gbm_seq_cluster6[1:5,1:5]

table(row.names(gbm_seq_cluster2)==row.names(gbm_seq_cluster5))
table(row.names(gbm_seq_cluster5)==row.names(gbm_seq_cluster6))

caf_like_cells <- cbind(gbm_seq_cluster2, gbm_seq_cluster5, gbm_seq_cluster6)
caf_like_cells_t <- t(caf_like_cells)
caf_like_cells_t_wanted <- subset(caf_like_cells_t, caf_like_cells_t[,'ACTA2'] > 0 & caf_like_cells_t[,'PTPRC'] == 0 & caf_like_cells_t[,'RGS5'] == 0 & caf_like_cells_t[,'GFAP'] == 0)
dim(caf_like_cells_t_wanted)
row.names(caf_like_cells_t_wanted)
interest <- sample(row.names(cluster_want_one), 119)
for(i in interest){
  sc_umap@cluster$kpart[names(sc_umap@cluster$kpart) == i] = 174
}

interest <- sample(row.names(cluster_want_two), 60)
for(i in interest){
  sc_umap@cluster$kpart[names(sc_umap@cluster$kpart) == i] = 175
}

interest <- sample(row.names(caf_like_cells_t_wanted), 278)
for(i in interest){
  sc_umap@cluster$kpart[names(sc_umap@cluster$kpart) == i] = 172
}

interest <- sample(row.names(vascular_cells), 115)
for(i in interest){
  sc_umap@cluster$kpart[names(sc_umap@cluster$kpart) == i] = 173
}

test <- sc_umap@cluster$kpart
test <- data.frame(test)
test$test[test$test ==1] <- "Cluster 1"
test$test[test$test ==2] <- "Cluster 2"
test$test[test$test ==3] <- "Cluster 3"
test$test[test$test ==4] <- "Cluster 4"
test$test[test$test ==5] <- "Cluster 5"
test$test[test$test ==6] <- "Cluster 6"
test$test[test$test ==7] <- "Cluster 7"
test$test[test$test ==8] <- "Cluster 8"
test$test[test$test ==173] <- "Cluster 9"
test$test[test$test ==174] <- "Cluster 10"
test$test[test$test ==175] <- "Cluster 11"
test$test[test$test ==172] <- "Cluster 12"
test$test <- factor(test$test, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "Cluster 11", "Cluster 12"))

gbm_seq_full <- read.table('/home/phillip/All_normal_cells_33K_GBM_seq.txt', header = T, row.names = 1)
dim(gbm_seq_full)
gbm_seq_full <- gbm_seq_full[,order(colnames(gbm_seq_full))]
gbm_seq_full <- t(gbm_seq_full)

meta_data <- read.csv(file.choose(), header = T, row.names = 1) #read:meta.csv
meta_data_normal <- subset(meta_data, meta_data$Tumor.Normal.Classification == 'Normal'| meta_data$Tumor.Normal.Classification == '')
meta_data_normal <- meta_data_normal[order(row.names(meta_data_normal)),]

gbm_seq_full_final <- gbm_seq_full[match(row.names(test), row.names(gbm_seq_full)),]
table(sort(row.names(gbm_seq_full_final) == row.names(test)))
meta_data_normal_final <- meta_data_normal[match(row.names(gbm_seq_full_final), row.names(meta_data_normal)),]
table(sort(row.names(gbm_seq_full_final) == row.names(meta_data_normal_final)))
table(sort(row.names(meta_data_normal_final) == row.names(test)))
wanted <- c('FN1', 'ACTA2', 'TAGLN')
gbm_seq_full_final <- gbm_seq_full_final[,match(wanted, colnames(gbm_seq_full_final))]
final <- cbind(row.names(test), test, row.names(gbm_seq_full_final), gbm_seq_full_final, row.names(meta_data_normal_final), meta_data_normal_final)
write.table(final, file='/home/phillip/matrix_for_bubble_pot.txt', col.names = T, row.names = T, quote = T)

matrix_for_bubble_plot <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/FN1_Miscellaneous_Figures/scRNAseq_BubblePlot/matrix_for_bubble_plot_to_analyze.txt', header = T, row.names = 1, sep = '\t') #read:matrix_for_bubble_plot_to_analyze
dim(matrix_for_bubble_plot)
sum(matrix_for_bubble_plot$FN1)

matrix_for_bubble_plot_caf <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'CAF')
(sum(matrix_for_bubble_plot_caf$FN1)/2312415)*100 #7.11%
mean(matrix_for_bubble_plot_caf$FN1) #879.359

matrix_for_bubble_plot_Bcell <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'B Cells')
(sum(matrix_for_bubble_plot_Bcell$FN1)/2312415)*100 #0.14%
mean(matrix_for_bubble_plot_Bcell$FN1) #2.356546

matrix_for_bubble_plot_CGEiN <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'CGE iN')
(sum(matrix_for_bubble_plot_CGEiN$FN1)/2312415)*100 #0.009%
mean(matrix_for_bubble_plot_CGEiN$FN1) #6.002569

matrix_for_bubble_plot_DividingBCells <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Dividing B Cells')
(sum(matrix_for_bubble_plot_DividingBCells$FN1)/2312415)*100 #0.07%
mean(matrix_for_bubble_plot_DividingBCells$FN1) #3.561877

matrix_for_bubble_plot_DividingNeuron <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Dividing Neuron')
(sum(matrix_for_bubble_plot_DividingNeuron$FN1)/2312415)*100 #1.4%
mean(matrix_for_bubble_plot_DividingNeuron$FN1) #38.804441

matrix_for_bubble_plot_DividingOPC <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Dividing OPC')
(sum(matrix_for_bubble_plot_DividingOPC$FN1)/2312415)*100 #0.03%
mean(matrix_for_bubble_plot_DividingOPC$FN1) #152.9186

matrix_for_bubble_plot_DividingProgenitor <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Dividing Progenitor')
(sum(matrix_for_bubble_plot_DividingProgenitor$FN1)/2312415)*100 #1.9%
mean(matrix_for_bubble_plot_DividingProgenitor$FN1) #178.184

matrix_for_bubble_plot_Endothelial<- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Endothelial')
(sum(matrix_for_bubble_plot_Endothelial$FN1)/2312415)*100 #23.3%
mean(matrix_for_bubble_plot_Endothelial$FN1) #675.8977

matrix_for_bubble_plot_GlycolyticProgenitor<- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Glycolytic Progenitor')
(sum(matrix_for_bubble_plot_GlycolyticProgenitor$FN1)/2312415)*100 #2.12%
mean(matrix_for_bubble_plot_GlycolyticProgenitor$FN1) #130.6346

matrix_for_bubble_plot_ImmatureAstrocyte <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Immature Astrocyte')
(sum(matrix_for_bubble_plot_ImmatureAstrocyte$FN1)/2312415)*100 #0.87%
mean(matrix_for_bubble_plot_ImmatureAstrocyte$FN1) #59.50015

matrix_for_bubble_plot_MatureIPCNewbornNeuron <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Mature IPC/Newborn Neuron')
(sum(matrix_for_bubble_plot_MatureIPCNewbornNeuron$FN1)/2312415)*100 #0.11%
mean(matrix_for_bubble_plot_MatureIPCNewbornNeuron$FN1) #162.2367

matrix_for_bubble_plot_Microglia <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Microglia')
(sum(matrix_for_bubble_plot_Microglia$FN1)/2312415)*100 #6.9%
mean(matrix_for_bubble_plot_Microglia$FN1) #87.94581

matrix_for_bubble_plot_MixedProgenitorNeuron <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Mixed Progenitor/Neuron')
(sum(matrix_for_bubble_plot_MixedProgenitorNeuron$FN1)/2312415)*100 #1.29%
mean(matrix_for_bubble_plot_MixedProgenitorNeuron$FN1) #102.8793

matrix_for_bubble_plot_Neuron <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Neuron')
(sum(matrix_for_bubble_plot_Neuron$FN1)/2312415)*100 #0.29%
mean(matrix_for_bubble_plot_Neuron$FN1) #16.75783

matrix_for_bubble_plot_Oligodendroycte <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Oligodendroycte')
(sum(matrix_for_bubble_plot_Oligodendroycte$FN1)/2312415)*100 #0.0%
mean(matrix_for_bubble_plot_Oligodendroycte$FN1) #0

matrix_for_bubble_plot_OPC <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'OPC')
(sum(matrix_for_bubble_plot_OPC$FN1)/2312415)*100 #0.0%
mean(matrix_for_bubble_plot_OPC$FN1) #0

matrix_for_bubble_plot_Pericyte <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Pericyte')
(sum(matrix_for_bubble_plot_Pericyte$FN1)/2312415)*100 #2.0%
mean(matrix_for_bubble_plot_Pericyte$FN1) #819.9982

matrix_for_bubble_plot_ProteoplasmicAstrocyte <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Protoplasmic Astrocyte')
(sum(matrix_for_bubble_plot_ProteoplasmicAstrocyte$FN1)/2312415)*100 #0.4%
mean(matrix_for_bubble_plot_ProteoplasmicAstrocyte$FN1) #233.4528

matrix_for_bubble_plot_RadialGlia <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Radial Glia')
(sum(matrix_for_bubble_plot_RadialGlia$FN1)/2312415)*100 #50.3%
mean(matrix_for_bubble_plot_RadialGlia$FN1) #536.2353

matrix_for_bubble_plot_RedBloodCells <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Red blood cells')
(sum(matrix_for_bubble_plot_RedBloodCells$FN1)/2312415)*100 #0.0%
mean(matrix_for_bubble_plot_RedBloodCells$FN1) #0

matrix_for_bubble_plot_TumorAssociatedMacrophage <- subset(matrix_for_bubble_plot, matrix_for_bubble_plot$Cell.Type.Assignment.RaceID3 == 'Tumor Associated Macrophage')
(sum(matrix_for_bubble_plot_TumorAssociatedMacrophage$FN1)/2312415)*100 #0.9%
mean(matrix_for_bubble_plot_TumorAssociatedMacrophage$FN1) #224.3326

fn1_df <- data.frame(Cell_Type = c('CAF',	'Pericyte',	'Endothelial',	'Radial Glia',	'Protoplasmic Astrocyte',	'Tumor Associated Macrophage',	'Dividing Progenitor',	'Mature IPC/Newborn Neuron',	'Dividing OPC',	'Glycolytic Progenitor',	'Mixed Progenitor/Neuron',	'Microglia',	'Immature Astrocyte',	'Dividing Neuron',	'Neuron',	'CGE iN',	'Dividing B Cells',	'B Cells',	'Oligodendrocyte',	'OPC',	'Red Blood Cells'), Percent_of_Total_Expression = c(7.11,	2,	23.3,	50.3,	0.4,	0.9,	1.9,	0.11,	0.03,	2.12,	1.29,	6.9,	0.87,	1.4,	0.29,	0.009,	0.07,	0.14,	0, 0,	0), Average_Gene_Expression = c(879.359,	819.9982,	675.8977,	536.2353,	233.4528,	224.3326,	178.184,	162.2367,	152.9186,	130.6346,	102.8793,	87.94581,	59.50015,	38.804441,	16.75783,	6.002569,	3.561877,	2.356546,	0,	0,	0))
fn1_df <- fn1_df[order(fn1_df$Average_Gene_Expression),]
fn1_df$Gene_Expression <- fn1_df$Average_Gene_Expression
mid <- mean(fn1_df$Gene_Expression)
ggplot(fn1_df, aes(x=reorder(Cell_Type, -Average_Gene_Expression), y=Average_Gene_Expression, color = Gene_Expression, size = Percent_of_Total_Expression)) +
  geom_point(alpha=5) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "grey", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 20), axis.text.y = element_text(face = "bold", size = 20), axis.title.x = element_text(size=22, face="bold", colour = "black"), axis.title.y = element_text(size=22, face="bold", colour = "black"), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black")) + 
  labs(x = "Cell Type", y = "FN1 Gene Expression", size = "Percent Expression", color = "Gene Expression") 

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=5) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 




