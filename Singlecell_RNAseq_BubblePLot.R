library(pheatmap)
library(limma)
library(ggplot2)
library(ggrepel)
library(plotly)
library(EnhancedVolcano)
library(dplyr)
library(tidyr)
library(devtools)
require(ggpubr)

#Marker Gene Bubble Plot
#Single Cell RNA-seq bubble plot (initial data matrix genreated on Deyou's computer)
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
final <- cbind(row.names(test), test, row.names(meta_data_normal_final), meta_data_normal_final, row.names(gbm_seq_full_final), gbm_seq_full_final)
write.table(final, file='/home/phillip/matrix_for_bubble_pot_figure1a.txt', col.names = T, row.names = T, quote = T)

#Make Bubble Plot in R
library(pheatmap)
table <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/ID_CAF_scRNAseq/matrix_for_bubble_pot_figure1a.txt', header = T, row.names = 1, sep = ' ')
dim(table)
table_test <- table
table_test$Cell.Type.Assignment[table_test$test == 'Cluster 12'] <- 'CAF'
table(sort(table_test$Cell.Type.Assignment))
genes_wanted <- c('ACTA2', 'VIM', 'LOX', 'CAV1', 'COL1A1', 'COL4A1', 'COL5A1', 'COL6A1', 'PTPRC', 'CD79A', 'CD79B', 'GAD2', 'ISL1', 'CCL22', 'FBXO43', 'NR2F2', 'SPC25', 'NMU', 'PENK', 'DEPDC1', 'CKAP2L', 'PECAM1', 'VWF', 'CDH5', 'NDUFA4L2', 'NUPR1', 'FXYD5', 'GFAP', 'ALDH1L1', 'ADCYAP1', 'BHLHE22', 'CSF1R', 'C1QC', 'PODXL', 'SLN', 'CNTNAP2', 'ELAVL4', 'PCSK6', 'TPPP', 'C1QTNF4', 'ITM2A', 'RGS5', 'PRELP', 'SLC7A10', 'GABRG1', 'CHI3L1', 'S100A16', 'HBD', 'HBM', 'CD14', 'CXCL1')
table_genes_wanted <- table_test[,match(genes_wanted, colnames(table_test))]
table_genes_wanted <- na.omit(table_genes_wanted)
dim(table_genes_wanted)

table_final <- table_test[match(row.names(table_genes_wanted), row.names(table_test)),]
table(sort(row.names(table_final) == row.names(table_genes_wanted)))
table(sort(table_final$Cell.Type.Assignment))
table_genes_wanted_final <- cbind(table_final$Cell.Type.Assignment, table_genes_wanted)

my_sample_col <- data.frame(cell_type = rep(c(table_final$Cell.Type.Assignment), 1))
row.names(my_sample_col) <- row.names(table_final)
class(my_sample_col)
pheatmap(table_genes_wanted, scale = "none", cluster_cols = F,  show_rownames = F, annotation_row = my_sample_col)

matrix_needed <- table_genes_wanted_final %>%
  group_by(table_genes_wanted_final$`table_final$Cell.Type.Assignment`) %>%
  summarise_if(is.numeric, mean)

matrix_needed <- data.frame(matrix_needed)  

test <- matrix_needed %>%
  select(table_genes_wanted_final..table_final.Cell.Type.Assignment.:CXCL1) %>%
  gather(Gene, Average_Expression, -table_genes_wanted_final..table_final.Cell.Type.Assignment.)

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = Average_Expression,color)) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Average Expression") 

#Obtain percentage of cells that express the marker genes
table_genes_wanted_final_caf <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'CAF')
dim(table_genes_wanted_final_caf)
write.table((colSums(table_genes_wanted_final_caf != 0)/187)*100)

table_genes_wanted_final_Bcells <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'B Cells')
dim(table_genes_wanted_final_Bcells)
write.table((colSums(table_genes_wanted_final_Bcells != 0)/1334)*100)

table_genes_wanted_final_CGEin <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'CGE iN')
dim(table_genes_wanted_final_CGEin)
write.table((colSums(table_genes_wanted_final_CGEin != 0)/35)*100)

table_genes_wanted_final_Dividing_B_cells <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing B Cells')
dim(table_genes_wanted_final_Dividing_B_cells)
write.table((colSums(table_genes_wanted_final_Dividing_B_cells != 0)/455)*100)

table_genes_wanted_final_Dividing_Neurons <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing Neuron')
dim(table_genes_wanted_final_Dividing_Neurons)
write.table((colSums(table_genes_wanted_final_Dividing_B_cells != 0)/845)*100)

table_genes_wanted_final_Dividing_Progenitor <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing OPC')
dim(table_genes_wanted_final_Dividing_Progenitor)
write.table((colSums(table_genes_wanted_final_Dividing_Progenitor != 0)/250)*100)

table_genes_wanted_final_Dividing_Progenitor <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing Progenitor')
dim(table_genes_wanted_final_Dividing_Progenitor)
write.table((colSums(table_genes_wanted_final_Dividing_Progenitor != 0)/250)*100)

table_genes_wanted_final_Endothelial <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Endothelial')
dim(table_genes_wanted_final_Endothelial)
write.table((colSums(table_genes_wanted_final_Endothelial != 0)/798)*100)

table_genes_wanted_final_Glycolytic_Progenitor <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Glycolytic Progenitor')
dim(table_genes_wanted_final_Glycolytic_Progenitor)
write.table((colSums(table_genes_wanted_final_Glycolytic_Progenitor != 0)/375)*100)

table_genes_wanted_final_Immature_astrocyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Immature Astrocyte')
dim(table_genes_wanted_final_Immature_astrocyte)
write.table((colSums(table_genes_wanted_final_Immature_astrocyte != 0)/340)*100)

table_genes_wanted_final_MatureIPC <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Mature IPC/Newborn Neuron')
dim(table_genes_wanted_final_MatureIPC)
write.table((colSums(table_genes_wanted_final_MatureIPC != 0)/15)*100)

table_genes_wanted_final_Microglia <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Microglia')
dim(table_genes_wanted_final_Microglia)
write.table((colSums(table_genes_wanted_final_Microglia != 0)/1829)*100)

table_genes_wanted_final_MixedProgenitorNeuron <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Mixed Progenitor/Neuron')
dim(table_genes_wanted_final_MixedProgenitorNeuron)
write.table((colSums(table_genes_wanted_final_MixedProgenitorNeuron != 0)/289)*100)

table_genes_wanted_final_Neuron <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Neuron')
dim(table_genes_wanted_final_Neuron)
write.table((colSums(table_genes_wanted_final_Neuron != 0)/405)*100)

table_genes_wanted_final_Oligodendrocyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Oligodendrocyte')
dim(table_genes_wanted_final_Oligodendrocyte)
write.table((colSums(table_genes_wanted_final_Oligodendrocyte != 0)/337)*100)

table_genes_wanted_final_OPC <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'OPC')
dim(table_genes_wanted_final_OPC)
write.table((colSums(table_genes_wanted_final_OPC != 0)/17)*100)

table_genes_wanted_final_Pericyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Pericyte')
dim(table_genes_wanted_final_Pericyte)
write.table((colSums(table_genes_wanted_final_Pericyte != 0)/56)*100)

table_genes_wanted_final_ProtoplasmicAstrocyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Protoplasmic Astrocyte')
dim(table_genes_wanted_final_ProtoplasmicAstrocyte)
write.table((colSums(table_genes_wanted_final_ProtoplasmicAstrocyte != 0)/42)*100)

table_genes_wanted_final_RadialGlia <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Radial Glia')
dim(table_genes_wanted_final_RadialGlia)
write.table((colSums(table_genes_wanted_final_RadialGlia != 0)/2167)*100)

table_genes_wanted_final_RedBloodCells <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Red blood cells')
dim(table_genes_wanted_final_RedBloodCells)
write.table((colSums(table_genes_wanted_final_RedBloodCells != 0)/24)*100)

table_genes_wanted_final_TumorAssociatedMacropahge <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Tumor Associated Macrophage')
dim(table_genes_wanted_final_TumorAssociatedMacropahge)
write.table((colSums(table_genes_wanted_final_TumorAssociatedMacropahge != 0)/93)*100)

single_cell_percent_expression_dataset <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/ID_CAF_scRNAseq/SingleCell_Average_Expressing_Cells.txt', header = T, sep = '\t')

test.2 <- single_cell_percent_expression_dataset %>%
  select(Cell.Type.Assignment:CXCL1) %>%
  gather(Gene, PercentExpression, -Cell.Type.Assignment)
head(test)
head(test.2)
table(sort(test$table_genes_wanted_final..table_final.Cell.Type.Assignment.)==test.2$Cell.Type.Assignment)

test$PercentExpression <- test.2$PercentExpression 

mid <- mean(log2(test$Average_Expression+1))
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient(low="grey", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), legend.title = element_text(size=12, face="bold", colour = "black"), legend.text = element_text(size=12, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Expression+1)") 

mid <- mean(test$Average_Expression+1)
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 


##################################################################################
#Run script below to generate bubble plot
library(pheatmap)
table <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/ID_CAF_scRNAseq/matrix_for_bubble_pot_figure1a.txt', header = T, row.names = 1, sep = ' ')
dim(table)
table_test <- table
table_test$Cell.Type.Assignment[table_test$test == 'Cluster 12'] <- 'CAF'
table(sort(table_test$Cell.Type.Assignment))
genes_wanted <- c('ACTA2', 'VIM', 'LOX', 'CAV1', 'COL1A1', 'COL4A1', 'COL5A1', 'COL6A1', 'PTPRC', 'CD79A', 'CD79B', 'GAD2', 'ISL1', 'CCL22', 'FBXO43', 'NR2F2', 'SPC25', 'NMU', 'PENK', 'DEPDC1', 'CKAP2L', 'PECAM1', 'VWF', 'CDH5', 'NDUFA4L2', 'NUPR1', 'FXYD5', 'GFAP', 'ALDH1L1', 'ADCYAP1', 'BHLHE22', 'CSF1R', 'C1QC', 'PODXL', 'SLN', 'CNTNAP2', 'ELAVL4', 'PCSK6', 'TPPP', 'C1QTNF4', 'ITM2A', 'RGS5', 'PRELP', 'SLC7A10', 'GABRG1', 'CHI3L1', 'S100A16', 'HBD', 'HBM', 'CD14', 'CXCL1')
table_genes_wanted <- table_test[,match(genes_wanted, colnames(table_test))]
table_genes_wanted <- na.omit(table_genes_wanted)
dim(table_genes_wanted)

table_final <- table_test[match(row.names(table_genes_wanted), row.names(table_test)),]
table(sort(row.names(table_final) == row.names(table_genes_wanted)))
table(sort(table_final$Cell.Type.Assignment))
table_genes_wanted_final <- cbind(table_final$Cell.Type.Assignment, table_genes_wanted)

matrix_needed <- table_genes_wanted_final %>%
  group_by(table_genes_wanted_final$`table_final$Cell.Type.Assignment`) %>%
  summarise_if(is.numeric, mean)

matrix_needed <- data.frame(matrix_needed)  

test <- matrix_needed %>%
  select(table_genes_wanted_final..table_final.Cell.Type.Assignment.:CXCL1) %>%
  gather(Gene, Average_Expression, -table_genes_wanted_final..table_final.Cell.Type.Assignment.)

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = Average_Expression,color)) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Average Expression") 

#Obtain percentage of cells that express the marker genes
table_genes_wanted_final_caf <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'CAF')
dim(table_genes_wanted_final_caf)
write.table((colSums(table_genes_wanted_final_caf != 0)/187)*100)

table_genes_wanted_final_Bcells <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'B Cells')
dim(table_genes_wanted_final_Bcells)
write.table((colSums(table_genes_wanted_final_Bcells != 0)/1334)*100)

table_genes_wanted_final_CGEin <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'CGE iN')
dim(table_genes_wanted_final_CGEin)
write.table((colSums(table_genes_wanted_final_CGEin != 0)/35)*100)

table_genes_wanted_final_Dividing_B_cells <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing B Cells')
dim(table_genes_wanted_final_Dividing_B_cells)
write.table((colSums(table_genes_wanted_final_Dividing_B_cells != 0)/455)*100)

table_genes_wanted_final_Dividing_Neurons <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing Neuron')
dim(table_genes_wanted_final_Dividing_Neurons)
write.table((colSums(table_genes_wanted_final_Dividing_B_cells != 0)/845)*100)

table_genes_wanted_final_Dividing_Progenitor <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing OPC')
dim(table_genes_wanted_final_Dividing_Progenitor)
write.table((colSums(table_genes_wanted_final_Dividing_Progenitor != 0)/250)*100)

table_genes_wanted_final_Dividing_Progenitor <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Dividing Progenitor')
dim(table_genes_wanted_final_Dividing_Progenitor)
write.table((colSums(table_genes_wanted_final_Dividing_Progenitor != 0)/250)*100)

table_genes_wanted_final_Endothelial <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Endothelial')
dim(table_genes_wanted_final_Endothelial)
write.table((colSums(table_genes_wanted_final_Endothelial != 0)/798)*100)

table_genes_wanted_final_Glycolytic_Progenitor <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Glycolytic Progenitor')
dim(table_genes_wanted_final_Glycolytic_Progenitor)
write.table((colSums(table_genes_wanted_final_Glycolytic_Progenitor != 0)/375)*100)

table_genes_wanted_final_Immature_astrocyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Immature Astrocyte')
dim(table_genes_wanted_final_Immature_astrocyte)
write.table((colSums(table_genes_wanted_final_Immature_astrocyte != 0)/340)*100)

table_genes_wanted_final_MatureIPC <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Mature IPC/Newborn Neuron')
dim(table_genes_wanted_final_MatureIPC)
write.table((colSums(table_genes_wanted_final_MatureIPC != 0)/15)*100)

table_genes_wanted_final_Microglia <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Microglia')
dim(table_genes_wanted_final_Microglia)
write.table((colSums(table_genes_wanted_final_Microglia != 0)/1829)*100)

table_genes_wanted_final_MixedProgenitorNeuron <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Mixed Progenitor/Neuron')
dim(table_genes_wanted_final_MixedProgenitorNeuron)
write.table((colSums(table_genes_wanted_final_MixedProgenitorNeuron != 0)/289)*100)

table_genes_wanted_final_Neuron <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Neuron')
dim(table_genes_wanted_final_Neuron)
write.table((colSums(table_genes_wanted_final_Neuron != 0)/405)*100)

table_genes_wanted_final_Oligodendrocyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Oligodendrocyte')
dim(table_genes_wanted_final_Oligodendrocyte)
write.table((colSums(table_genes_wanted_final_Oligodendrocyte != 0)/337)*100)

table_genes_wanted_final_OPC <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'OPC')
dim(table_genes_wanted_final_OPC)
write.table((colSums(table_genes_wanted_final_OPC != 0)/17)*100)

table_genes_wanted_final_Pericyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Pericyte')
dim(table_genes_wanted_final_Pericyte)
write.table((colSums(table_genes_wanted_final_Pericyte != 0)/56)*100)

table_genes_wanted_final_ProtoplasmicAstrocyte <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Protoplasmic Astrocyte')
dim(table_genes_wanted_final_ProtoplasmicAstrocyte)
write.table((colSums(table_genes_wanted_final_ProtoplasmicAstrocyte != 0)/42)*100)

table_genes_wanted_final_RadialGlia <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Radial Glia')
dim(table_genes_wanted_final_RadialGlia)
write.table((colSums(table_genes_wanted_final_RadialGlia != 0)/2167)*100)

table_genes_wanted_final_RedBloodCells <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Red blood cells')
dim(table_genes_wanted_final_RedBloodCells)
write.table((colSums(table_genes_wanted_final_RedBloodCells != 0)/24)*100)

table_genes_wanted_final_TumorAssociatedMacropahge <- subset(table_genes_wanted_final, table_genes_wanted_final$`table_final$Cell.Type.Assignment` == 'Tumor Associated Macrophage')
dim(table_genes_wanted_final_TumorAssociatedMacropahge)
write.table((colSums(table_genes_wanted_final_TumorAssociatedMacropahge != 0)/93)*100)

single_cell_percent_expression_dataset <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/ID_CAF_scRNAseq/SingleCell_Percent_Expressing_Cells.txt', header = T, sep = '\t')

test.2 <- single_cell_percent_expression_dataset %>%
  select(Cell.Type.Assignment:CXCL1) %>%
  gather(Gene, PercentExpression, -Cell.Type.Assignment)
head(test)
head(test.2)
table(sort(test$table_genes_wanted_final..table_final.Cell.Type.Assignment.)==test.2$Cell.Type.Assignment)

test$PercentExpression <- test.2$PercentExpression 

mid <- mean(log2(test$Average_Expression+1))
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient(low="grey", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), legend.title = element_text(size=12, face="bold", colour = "black"), legend.text = element_text(size=12, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Expression+1)") 

mid <- mean(test$Average_Expression+1)
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 

mid <- mean(test$Average_Expression+1)
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = ifelse(PercentExpression==0, NA,PercentExpression), color =(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
#Run the below script
#Ultimate script to make bubble plot
library(pheatmap)
table <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/ID_CAF_scRNAseq/matrix_for_bubble_pot_figure1a.txt', header = T, row.names = 1, sep = ' ')
dim(table)
table_test <- table
table_test$Cell.Type.Assignment[table_test$test == 'Cluster 12'] <- 'CAF'
table(sort(table_test$Cell.Type.Assignment))
genes_wanted <- c('ACTA2', 'VIM', 'LOX', 'CAV1', 'COL1A1', 'COL4A1', 'COL5A1', 'COL6A1', 'PTPRC', 'CD79A', 'CD79B', 'GAD2', 'ISL1', 'CCL22', 'FBXO43', 'NR2F2', 'SPC25', 'NMU', 'PENK', 'DEPDC1', 'CKAP2L', 'PECAM1', 'VWF', 'CDH5', 'NDUFA4L2', 'NUPR1', 'FXYD5', 'GFAP', 'ALDH1L1', 'ADCYAP1', 'BHLHE22', 'CSF1R', 'C1QC', 'PODXL', 'SLN', 'CNTNAP2', 'ELAVL4', 'PCSK6', 'TPPP', 'C1QTNF4', 'ITM2A', 'RGS5', 'PRELP', 'SLC7A10', 'GABRG1', 'CHI3L1', 'S100A16', 'HBD', 'HBM', 'CD14', 'CXCL1')
table_genes_wanted <- table_test[,match(genes_wanted, colnames(table_test))]
table_genes_wanted <- na.omit(table_genes_wanted)
dim(table_genes_wanted)

table_final <- table_test[match(row.names(table_genes_wanted), row.names(table_test)),]
table(sort(row.names(table_final) == row.names(table_genes_wanted)))
table(sort(table_final$Cell.Type.Assignment))
table_genes_wanted_final <- cbind(table_final$Cell.Type.Assignment, table_genes_wanted)

my_sample_col <- data.frame(cell_type = rep(c(table_final$Cell.Type.Assignment), 1))
row.names(my_sample_col) <- row.names(table_final)
class(my_sample_col)
pheatmap(table_genes_wanted, scale = "none", cluster_cols = F,  show_rownames = F, annotation_row = my_sample_col)

matrix_needed <- table_genes_wanted_final %>%
  group_by(table_genes_wanted_final$`table_final$Cell.Type.Assignment`) %>%
  summarise_if(is.numeric, mean)

matrix_needed <- data.frame(matrix_needed)  

test <- matrix_needed %>%
  select(table_genes_wanted_final..table_final.Cell.Type.Assignment.:CXCL1) %>%
  gather(Gene, Average_Expression, -table_genes_wanted_final..table_final.Cell.Type.Assignment.)

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = Average_Expression,color)) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Average Expression") 

single_cell_percent_expression_dataset <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/ID_CAF_scRNAseq/SingleCell_Percent_Expressing_Cells.txt', header = T, sep = '\t')

test.2 <- single_cell_percent_expression_dataset %>%
  select(Cell.Type.Assignment:CXCL1) %>%
  gather(Gene, PercentExpression, -Cell.Type.Assignment)
head(test)
head(test.2)
table(sort(test$table_genes_wanted_final..table_final.Cell.Type.Assignment.)==test.2$Cell.Type.Assignment)

test$PercentExpression <- test.2$PercentExpression 

mid <- mean(log2(test$Average_Expression+1))
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 

ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1', 'PTPRC',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =log2(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient(low="grey", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), legend.title = element_text(size=12, face="bold", colour = "black"), legend.text = element_text(size=12, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Expression+1)") 

mid <- mean(test$Average_Expression+1)
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = PercentExpression, color =(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 

mid <- mean(test$Average_Expression+1)
ggplot(data.frame(test), aes(x= factor(Gene, level = c('ACTA2',	'VIM',	'LOX',	'CAV1',	'COL1A1',	'COL4A1',	'COL5A1',	'COL6A1',	'CD79A',	'CD79B',	'GAD2',	'ISL1',	'CCL22',	'FBXO43',	'NR2F2',	'SPC25',	'NMU',	'PENK',	'DEPDC1',	'CKAP2L',	'PECAM1', 'VWF', 'CDH5',	'NDUFA4L2',	'NUPR1',	'FXYD5',	'GFAP',	'ALDH1L1',	'ADCYAP1',	'BHLHE22',	'CSF1R',	'C1QC',	'PODXL',	'SLN',	'CNTNAP2',	'ELAVL4',	'PCSK6',	'TPPP',	'C1QTNF4',	'ITM2A',	'RGS5',	'PRELP',	'SLC7A10',	'GABRG1',	'CHI3L1',	'S100A16',	'HBD',	'HBM',	'CD14',	'CXCL1')), y= factor(table_genes_wanted_final..table_final.Cell.Type.Assignment., level = c('Tumor Associated Macrophage',	'Red blood cells',	'Radial Glia',	'Protoplasmic Astrocyte',	'Pericyte',	'OPC',	'Oligodendrocyte',	'Neuron',	'Mixed Progenitor/Neuron',		'Microglia',	'Mature IPC/Newborn Neuron',		'Immature Astrocyte',	'Glycolytic Progenitor',	'Endothelial',	'Dividing Progenitor',	'Dividing OPC',	'Dividing Neuron',	'Dividing B Cells',	'CGE iN',	'B Cells',	'CAF')), size = ifelse(PercentExpression==0, NA,PercentExpression), color =(Average_Expression+1))) +
  geom_point(alpha=0.7) +
  theme_classic() + 
  scale_color_gradient2(midpoint=mid, low="blue", mid = "white", high="red") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black")) + 
  labs(x = "Marker Gene", y = "Cell Type", size = "Percent Expression", color = "log2(Count+1)") 


