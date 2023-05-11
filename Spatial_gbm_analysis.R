#GBM Spatial Analysis
library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(tidyverse)
library(celltalker)
library(devtools)

#Too load data make a folder on desktop this will contain the .h5 file and separate folder that contains the .png file and the .json file 
#original analysis
gbm_spatial_data <- Load10X_Spatial('/Users/phillipgalbo/Desktop/Spatial_Information', filename = "Parent_Visium_Human_Glioblastoma_filtered_feature_bc_matrix.h5")
plot1 <- VlnPlot(gbm_spatial_data, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(gbm_spatial_data, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)
gbm_spatial_data <- SCTransform(gbm_spatial_data, assay = "Spatial", verbose = FALSE)
SpatialFeaturePlot(gbm_spatial_data, features = c("HPCA", "TTR"))
p1 <- SpatialFeaturePlot(gbm_spatial_data, features = "TTR", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(gbm_spatial_data, features = "TTR", alpha = c(0.1, 1))
p1 + p2
gbm_spatial_data <- RunPCA(gbm_spatial_data, assay = "SCT", verbose = FALSE)
gbm_spatial_data <- FindNeighbors(gbm_spatial_data, reduction = "pca", dims = 1:30)
#resolution 0.6 THIS is the best result
gbm_spatial_data_0.6 <- FindClusters(gbm_spatial_data, verbose = FALSE, resolution = 0.6)
gbm_spatial_data_0.6 <- RunUMAP(gbm_spatial_data_0.6, reduction = "pca", dims = 1:30)
DimPlot(gbm_spatial_data_0.6, reduction = "umap", label = TRUE)
SpatialDimPlot(gbm_spatial_data_0.6, label = TRUE, label.size = 3)
gbm_spatial_data_markers_0.6 <- FindAllMarkers(gbm_spatial_data_0.6, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
write.table(gbm_spatial_data_markers_0.6, file='/Users/phillipgalbo/Desktop/gbm_spatial_data_markers_0.6.txt', col.names = T, row.names = T, quote = F)
gbm_spatial_data_markers_0.6_top20 <- gbm_spatial_data_markers_0.6 %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
DoHeatmap(gbm_spatial_data_0.6, features = gbm_spatial_data_markers_0.6_top20$gene, size = 12) + NoLegend()
saveRDS(gbm_spatial_data_0.6, file = "/Users/phillipgalbo/Desktop/gbm_spatial_data_0.6.rds")

#CAF markers
gbm_spatial_data_0.6 <- readRDS('/Volumes/LaCie\ 1/GBM\ spatial\ dataset/resolution_0.6/Resolution\ 0.6/gbm_spatial_data_0.6.rds')
SpatialDimPlot(gbm_spatial_data_0.6, label = TRUE, label.size = 12) + NoLegend()
gbm_spatial_data_markers_0.6 <- FindAllMarkers(gbm_spatial_data_0.6, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
gbm_spatial_data_markers_0.6_top20 <- gbm_spatial_data_markers_0.6 %>% group_by(cluster) %>% top_n(n=20, wt = avg_log2FC)
DoHeatmap(gbm_spatial_data_0.6, features = gbm_spatial_data_markers_0.6_top20$gene) +  theme(legend.title = element_text(size=12, face="bold", colour = "black"), legend.text = element_text(size=12, face="bold", colour = "black"))

VlnPlot(gbm_spatial_data_0.6, feature = 'ACTA2')
VlnPlot(gbm_spatial_data_0.6, feature = 'PECAM1')
VlnPlot(gbm_spatial_data_0.6, feature = 'PECAM1')
VlnPlot(gbm_spatial_data_0.6, feature = 'RGS5')
VlnPlot(gbm_spatial_data_0.6, feature = 'PTPRC')


VlnPlot(gbm_spatial_data_0.6, feature = 'COL1A1')
VlnPlot(gbm_spatial_data_0.6, feature = 'COL1A2')
VlnPlot(gbm_spatial_data_0.6, feature = 'COL6A1')
VlnPlot(gbm_spatial_data_0.6, feature = 'TAGLN')
VlnPlot(gbm_spatial_data_0.6, feature = 'SLC39A14')

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "ACTA2") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "PECAM1") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "RGS5") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "PTPRC") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "GJB2") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "PTX3") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "SERPINE1") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "ACTG2") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "AKAP12") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "ITGA5") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "FN1") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
VlnPlot(gbm_spatial_data_0.6, feature = 'FN1') + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), legend.title = element_text(size=12, face="bold", colour = "black"), legend.text = element_text(size=12, face="bold", colour = "black"))

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "IL6") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CXCL8") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCL2") + theme(legend.key.size = unit(1, 'cm'), legend.title = element_text(size=20, face="bold", colour = "black"), legend.text = element_text(size=20, face="bold", colour = "black"))

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "COL1A1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "COL1A2")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "COL5A1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "COL6A1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "TAGLN")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CAV1")

#Endothelial/Pericyte
VlnPlot(gbm_spatial_data_0.6, feature = 'PECAM1')
VlnPlot(gbm_spatial_data_0.6, feature = 'RGS5')

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "PECAM1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "RGS5")

#Immune Cell
VlnPlot(gbm_spatial_data_0.6, feature = 'PTPRC')

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "PTPRC")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CD14")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CD79A")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CD79B")

#Cytokines and receptors
VlnPlot(gbm_spatial_data_0.6, feature = 'CXCL8')
VlnPlot(gbm_spatial_data_0.6, feature = 'IL6')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCL2')

SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CXCL8")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "IL6")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCL2")

#IL-8 Receptors
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CXCR1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CXCR2")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "ACKR1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "KDR")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "SDC1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "SDC2")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "SDC3")

VlnPlot(gbm_spatial_data_0.6, feature = 'CXCR1')
VlnPlot(gbm_spatial_data_0.6, feature = 'CXCR2')
VlnPlot(gbm_spatial_data_0.6, feature = 'ACKR1')
VlnPlot(gbm_spatial_data_0.6, feature = 'KDR')
VlnPlot(gbm_spatial_data_0.6, feature = 'SDC1')
VlnPlot(gbm_spatial_data_0.6, feature = 'SDC2')
VlnPlot(gbm_spatial_data_0.6, feature = 'SDC3')

#CCL2 Receptors
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "ACKR2")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "ACKR4")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCR1")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCR10")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCR2")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCR3")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCR4")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "CCR5")
SpatialFeaturePlot(gbm_spatial_data_0.6, features = "DARC")

VlnPlot(gbm_spatial_data_0.6, feature = 'ACKR2')
VlnPlot(gbm_spatial_data_0.6, feature = 'ACKR4')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCR1')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCR10')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCR2')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCR3')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCR4')
VlnPlot(gbm_spatial_data_0.6, feature = 'CCR5')
VlnPlot(gbm_spatial_data_0.6, feature = 'DARC')

#heatmap of different cell types; top five marker genes (cell type marker genes were derived from the outer radial glia-like cancer stem cells contirbute to heterogeneity of glioblastoma" manuscript)
marker_genes <- read.table(file.choose(), header = T, sep = '\t') #load file bhadurietal_marker_genes
marker_genes <- marker_genes[,1:7]

marker_genes_matureIPC_newborn_Neuron <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '1')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_matureIPC_newborn_Neuron$gene) + NoLegend()

marker_genes_glycolytic_progenitors <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '2' | marker_genes$cluster == '34')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_glycolytic_progenitors$gene) + NoLegend()

marker_genes_OPC <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '3' | marker_genes$cluster == '11' | marker_genes$cluster == '31' | marker_genes$cluster == '39')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_OPC$gene) + NoLegend()

marker_genes_oligodendrocyte <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '4' | marker_genes$cluster == '9' | marker_genes$cluster == '12' | marker_genes$cluster == '24' | marker_genes$cluster == '50')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_oligodendrocyte$gene) + NoLegend()

marker_genes_microglia <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '5' | marker_genes$cluster == '17' | marker_genes$cluster == '29' | marker_genes$cluster == '32' | marker_genes$cluster == '33' | marker_genes$cluster == '36' | marker_genes$cluster == '37' | marker_genes$cluster == '44' |marker_genes$cluster == '48')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_microglia$gene) + NoLegend()

marker_genes_mixed_progenitor_neuron <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '6' | marker_genes$cluster == '42')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_mixed_progenitor_neuron$gene) + NoLegend()

marker_genes_endothelial_cells <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '7' | marker_genes$cluster == '13' | marker_genes$cluster == '20' | marker_genes$cluster == '22')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_endothelial_cells$gene) + NoLegend()

marker_genes_protoplasmic_astrocyte <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '15')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_protoplasmic_astrocyte$gene) + NoLegend()

marker_genes_TAMs <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '8')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_TAMs$gene) + NoLegend()

marker_genes_pericyte <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '10' | marker_genes$cluster == '21')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_pericyte$gene) + NoLegend()

marker_genes_bcells <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '16')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_bcells$gene) + NoLegend()

marker_genes_radial_glia <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '19' | marker_genes$cluster == '23' | marker_genes$cluster == '26' | marker_genes$cluster == '27' | marker_genes$cluster == '28' | marker_genes$cluster == '38')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_radial_glia$gene) + NoLegend()

marker_genes_immature_astrocyte <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '41' | marker_genes$cluster == '43')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_immature_astrocyte$gene) + NoLegend()

marker_genes_neuron <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '45')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_neuron$gene) + NoLegend()

marker_genes_redbloodcell <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '46')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_redbloodcell$gene) + NoLegend()

marker_genes_CGEiN <- subset(marker_genes, marker_genes$p_val_adj < 0.05 & marker_genes$cluster == '47')
DoHeatmap(gbm_spatial_data_0.6, features = marker_genes_CGEiN$gene) + NoLegend()

pan_dcaf_markers <- c('MMP1','IGFBP3', 'MMP11', 'STC1', 'IGLL5', 'COL10A1', 'POSTN', 'CA12', 'THBS2', 'IDO1', 'CTHRC1', 'COMP', 'COL1A1', 'LRRC15', 'EVA1A', 'MFAP2', 'HAS1', 'SDC1', 'SLC38A5', 'LUM', 'COL11A1', 'TWIST2','COL3A1', 'SPON2', 'CXCL1', 'TREM1', 'TACSTD2', 'FN1', 'WNT2', 'SERPINE2', 'COL7A1', 'COL12A1', 'ANGPTL4', 'WISP1', 'COL1A2', 'CRABP2', 'SPOCK1', 'SERPINE1', 'DERL3', 'COL5A1', 'KCNN4', 'GBP5', 'LOXL2', 'TWIST1', 'PAPPA', 'SLAMF8', 'UCHL1', 'ADAM12', 'SYNDIG1', 'F2RL2', 'CTSK', 'PLAUR', 'HAS2', 'RSAD2', 'P4HA3', 'TNFAIP6', 'INHBA', 'RGCC', 'FNDC1', 'ASPN', 'STRA6')
DoHeatmap(gbm_spatial_data_0.6, features = pan_dcaf_markers) + NoLegend()

gbm_caf_markers <- read.table(file.choose(), sep = '\t') #read file titled: sc_CAF_gbm_markers.txt
DoHeatmap(gbm_spatial_data_0.6, features = gbm_caf_markers$V1) + NoLegend()

#Load Verhaak subtype gene signatures that are used for GSEA
verrhak_proneural <- read.table(file.choose(), header = T, sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = verrhak_proneural$VERHAAK_GLIOBLASTOMA_PRONEURAL) + NoLegend()
verrhak_classical <- read.table(file.choose(), header = T, sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = verrhak_classical$VERHAAK_GLIOBLASTOMA_CLASSICAL) + NoLegend()
verrhak_mesenchymal <- read.table(file.choose(), header = T, sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = verrhak_mesenchymal$VERHAAK_GLIOBLASTOMA_MESENCHYMAL) + NoLegend()

#Ivygap analysis
CT <- read.table(file.choose(), sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = CT$V1) + NoLegend()
CTmvp <- read.table(file.choose(), sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = CTmvp$V1) + NoLegend()
CTpan <- read.table(file.choose(), sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = CTpan$V1) + NoLegend()
LE <- read.table(file.choose(), sep = '\t')
DoHeatmap(gbm_spatial_data_0.6, features = LE$V1) + NoLegend()

#CellTalk analysis
gbm_spatial_data_0.6 <- readRDS(file.choose())

ramilowski_pairs$ligand <- gsub('IL8', 'CXCL8', ramilowski_pairs$ligand)
ramilowski_pairs$ligand <- gsub('IL6', 'CXCL6', ramilowski_pairs$ligand)

ligs <- as.character(unique(ramilowski_pairs$ligand))
recs <- as.character(unique(ramilowski_pairs$receptor))
ligs.present <- rownames(gbm_spatial_data_0.6)[rownames(gbm_spatial_data_0.6) %in% ligs]
recs.present <- rownames(gbm_spatial_data_0.6)[rownames(gbm_spatial_data_0.6) %in% recs]
genes.to.use <- union(ligs.present,recs.present)

markers <- FindAllMarkers(gbm_spatial_data_0.6, features=genes.to.use,only.pos=TRUE)
ligs.recs.use <- unique(markers$gene)

interactions.forward1 <- ramilowski_pairs[as.character(ramilowski_pairs$ligand) %in% ligs.recs.use,]
interactions.forward2 <- ramilowski_pairs[as.character(ramilowski_pairs$receptor) %in% ligs.recs.use,]
interact.for <- rbind(interactions.forward1,interactions.forward2)

expr.mat <- GetAssayData(gbm_spatial_data_0.6, slot="counts")
defined.clusters <- gbm_spatial_data_0.6@meta.data$seurat_clusters
defined.groups <- gbm_spatial_data_0.6@meta.data$seurat_clusters
defined.replicates <- gbm_spatial_data_0.6@meta.data$seurat_clusters

reshaped.matrices <- reshape_matrices(count.matrix=expr.mat, clusters=defined.clusters, groups=defined.groups, replicates=defined.replicates, ligands.and.receptors=interact.for)


consistent.lig.recs <- create_lig_rec_tib(exp.tib=reshaped.matrices,
                                          clusters=defined.clusters,
                                          groups=defined.groups,
                                          replicates=defined.replicates,
                                          cells.reqd=10,
                                          freq.pos.reqd=0.5,
                                          ligands.and.receptors=interact.for)

put.int<-putative_interactions(ligand.receptor.tibble=consistent.lig.recs,
                               clusters=defined.clusters,
                               groups=defined.groups,
                               freq.group.in.cluster=0.05,
                               ligands.and.receptors=interact.for)

#CellPhoneDB Dataset curation
gbm_spatial_data_0.6 <- readRDS(file.choose())
gbm_spatial_data_0.6_scale_data <- GetAssayData(object = gbm_spatial_data_0.6, slot = "scale.data")
gbm_spatial_data_0.6_count_data <- GetAssayData(object = gbm_spatial_data_0.6, slot = "counts")
write.table(gbm_spatial_data_0.6_scale_data, "/Users/phillipgalbo/Desktop/gbm_spatial_data_0.6_scale_data.txt", row.names=T, col.names = T, sep='\t')
meta_data <- data.frame(Cell=rownames(gbm_spatial_data_0.6@meta.data), cell_type=gbm_spatial_data_0.6@meta.data$seurat_clusters)
write.table(meta_data, "/Users/phillipgalbo/Desktop/meta_data.txt", row.names=F, col.names = T, sep='\t', quote = F)












