#GSVA of GBM subtypes between GBMs with high and low CAFs
library(genefilter)
library(GSVA)
library(Biobase)
library(reshape2)
library(stringr)
library(pheatmap)
library(DESeq2)
library(ggpubr)
library(patchwork)

gsva_gene_set <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/PMT_Preliminary_Data/GSVA_GBM_Subtype/GSVA_subtypes_matrix.csv', header = T)
list <- split(as.matrix(gsva_gene_set)[,1], gsva_gene_set[,2])

tcga_gbm_count <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/TCGA_GBM_counts.txt', header = T, row.names = 1, sep = '\t')
tcga_gbm_count <- 2^tcga_gbm_count
tcga_gbm_count[1:5,1:5]

tcga_gbm_gsva <- gsva(as.matrix(tcga_gbm_count), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T)
samples_wanted <- c('TCGA.02.0055.01A',	'TCGA.19.1389.02A',	'TCGA.14.0789.01A',	'TCGA.27.2524.01A',	'TCGA.28.2513.01A',	'TCGA.06.0190.01A',	'TCGA.14.0781.01B',	'TCGA.06.0645.01A',	'TCGA.06.0210.01A',	'TCGA.06.0139.01A',	'TCGA.06.0878.01A',	'TCGA.26.5139.01A',	'TCGA.06.5410.01A',	'TCGA.06.0644.01A',	'TCGA.06.0171.02A',	'TCGA.32.2615.01A',	'TCGA.32.4213.01A',	'TCGA.02.2486.01A',	'TCGA.06.2562.01A',	'TCGA.06.0168.01A',	'TCGA.41.4097.01A',	'TCGA.06.0211.01A',	'TCGA.41.3915.01A',	'TCGA.06.0156.01A',	'TCGA.06.5412.01A',	'TCGA.12.5299.01A',	'TCGA.27.1832.01A',	'TCGA.06.5418.01A',	'TCGA.76.4928.01B',	'TCGA.06.0187.01A',	'TCGA.06.0141.01A',	'TCGA.28.5213.01A',	'TCGA.28.5216.01A',	'TCGA.32.2638.01A',	'TCGA.19.2625.01A',	'TCGA.16.1045.01B',	'TCGA.12.0619.01A',	'TCGA.06.0750.01A',	'TCGA.06.0130.01A',	'TCGA.14.2554.01A',	'TCGA.16.0846.01A',	'TCGA.06.0157.01A',	'TCGA.02.2485.01A',	'TCGA.12.0616.01A',	'TCGA.06.2563.01A',	'TCGA.19.1787.01B',	'TCGA.06.2570.01A',	'TCGA.14.0871.01A',	'TCGA.28.5204.01A',	'TCGA.06.5417.01A',	'TCGA.26.5134.01A',	'TCGA.14.1825.01A',	'TCGA.06.0238.01A',	'TCGA.27.2521.01A',	'TCGA.19.2629.01A',	'TCGA.12.5295.01A',	'TCGA.41.5651.01A',	'TCGA.06.0744.01A',	'TCGA.32.1970.01A',	'TCGA.27.1835.01A',	'TCGA.12.3650.01A',	'TCGA.14.0790.01B',	'TCGA.28.5207.01A',	'TCGA.19.2624.01A',	'TCGA.28.1747.01C',	'TCGA.12.3653.01A',	'TCGA.19.1390.01A',	'TCGA.15.1444.01A',	'TCGA.26.1442.01A',	'TCGA.06.2557.01A',	'TCGA.06.0178.01A',	'TCGA.06.0749.01A',	'TCGA.06.2569.01A',	'TCGA.28.5220.01A',	'TCGA.26.5133.01A',	'TCGA.06.1804.01A',	'TCGA.06.5859.01A',	'TCGA.06.0221.02A',	'TCGA.19.5960.01A',	'TCGA.06.0129.01A')
tcga_gbm_gsva_wanted <- tcga_gbm_gsva[,match(samples_wanted, colnames(tcga_gbm_gsva))]
gene_sets_wanted <- c('VERHAAK_GLIOBLASTOMA_MESENCHYMAL', 'PHILIPS_MESENCHYMAL', 'NEFTEL_MES_META_MODULE', 'WANG_mGSC', 'VERHAAK_GLIOBLASTOMA_PRONEURAL', 'NEFTEL_NPC_META_MODULE', 'NEFTEL_OPC_META_MODULE', 'PHILIPS_PRONEURAL', 'WANG_pGSC', 'VERHAAK_GLIOBLASTOMA_CLASSICAL', 'NEFTEL_AC_META_MODULE')
tcga_gbm_gsva_wanted <- tcga_gbm_gsva_wanted[match(gene_sets_wanted, row.names(tcga_gbm_gsva_wanted)),]
my_sample_col <- data.frame(sample = rep(c("CAF_High", "CAF_Low"), c(40,40)))
row.names(my_sample_col) <- colnames(tcga_gbm_gsva_wanted)
my_colour = list(sample = c(CAF_High = "red", CAF_Low = "blue"))
pheatmap(tcga_gbm_gsva_wanted, scale = "row", cluster_cols = F, cluster_rows = F, color=colorRampPalette(c("blue","white","red"))(50), breaks=seq(-1, 1, length.out = 50),  show_rownames = T, show_colnames = F, annotation_col = my_sample_col, annotation_colors = my_colour)

cgga_325 <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/CGGA.325.gene.count.name.txt', header = T, row.names = 1, sep = '\t')
cgga_325_gbm_gsva <- gsva(as.matrix(cgga_325), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T)
samples_wanted <- c('CGGA_318',	'CGGA_1035',	'CGGA_859',	'CGGA_1177',	'CGGA_1045',	'CGGA_1342',	'CGGA_1240',	'CGGA_1124',	'CGGA_802',	'CGGA_624',	'CGGA_D26',	'CGGA_1338',	'CGGA_D36',	'CGGA_842',	'CGGA_D35',	'CGGA_1049',	'CGGA_1077',	'CGGA_1313',	'CGGA_604',	'CGGA_1393',	'CGGA_1073',	'CGGA_731',	'CGGA_D51',	'CGGA_1270',	'CGGA_1180',	'CGGA_D37',	'CGGA_876',	'CGGA_1015',	'CGGA_525',	'CGGA_1299',	'CGGA_837',	'CGGA_1083',	'CGGA_1224',	'CGGA_1105',	'CGGA_1114',	'CGGA_1001',	'CGGA_759',	'CGGA_1287',	'CGGA_1197',	'CGGA_1384',	'CGGA_1370',	'CGGA_D34',	'CGGA_D30',	'CGGA_899',	'CGGA_632',	'CGGA_1170',	'CGGA_1072',	'CGGA_808',	'CGGA_1039',	'CGGA_1129',	'CGGA_491',	'CGGA_1227',	'CGGA_J100',	'CGGA_1450',	'CGGA_494',	'CGGA_1409',	'CGGA_374',	'CGGA_1320',	'CGGA_804',	'CGGA_680',	'CGGA_1219',	'CGGA_1375',	'CGGA_1394',	'CGGA_669',	'CGGA_1251',	'CGGA_1283',	'CGGA_747',	'CGGA_1060',	'CGGA_1116',	'CGGA_822',	'CGGA_700',	'CGGA_518')
cgga_325_gbm_gsva_wanted <- cgga_325_gbm_gsva[,match(samples_wanted, colnames(cgga_325_gbm_gsva))]
gene_sets_wanted <- c('VERHAAK_GLIOBLASTOMA_MESENCHYMAL', 'PHILIPS_MESENCHYMAL', 'NEFTEL_MES_META_MODULE', 'WANG_mGSC', 'VERHAAK_GLIOBLASTOMA_PRONEURAL', 'NEFTEL_NPC_META_MODULE', 'NEFTEL_OPC_META_MODULE', 'PHILIPS_PRONEURAL', 'WANG_pGSC', 'VERHAAK_GLIOBLASTOMA_CLASSICAL', 'NEFTEL_AC_META_MODULE')
cgga_325_gbm_gsva_wanted <- cgga_325_gbm_gsva_wanted[match(gene_sets_wanted, row.names(cgga_325_gbm_gsva_wanted)),]
my_sample_col <- data.frame(sample = rep(c("CAF_High", "CAF_Low"), c(36,36)))
row.names(my_sample_col) <- colnames(cgga_325_gbm_gsva_wanted)
my_colour = list(sample = c(CAF_High = "red", CAF_Low = "blue"))
pheatmap(cgga_325_gbm_gsva_wanted, scale = "row", cluster_cols = F, cluster_rows = F, color=colorRampPalette(c("blue","white","red"))(50), breaks=seq(-1, 1, length.out = 50),  show_rownames = T, show_colnames = F, annotation_col = my_sample_col, annotation_colors = my_colour)

cgga_693 <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/CGGA.693.gene.count.name.txt', header = T, row.names = 1, sep = '\t')
cgga_693_gbm_gsva <- gsva(as.matrix(cgga_693), list, method = 'ssgsea', kcdf = 'Poisson', abs.ranking = T)
samples_wanted <- c('CGGA_P7',	'CGGA_1208',	'CGGA_P143',	'CGGA_1134',	'CGGA_1415',	'CGGA_1546',	'CGGA_1946',	'CGGA_1706',	'CGGA_1551',	'CGGA_1586',	'CGGA_869',	'CGGA_2088',	'CGGA_2075',	'CGGA_1972',	'CGGA_P182',	'CGGA_1106',	'CGGA_P89',	'CGGA_1744',	'CGGA_1612',	'CGGA_1492',	'CGGA_1480',	'CGGA_1451',	'CGGA_1503',	'CGGA_1814',	'CGGA_P610',	'CGGA_1441',	'CGGA_1708',	'CGGA_1901',	'CGGA_1626',	'CGGA_P585',	'CGGA_1687',	'CGGA_1538',	'CGGA_P499',	'CGGA_2078',	'CGGA_1452',	'CGGA_1481',	'CGGA_1505',	'CGGA_831',	'CGGA_1601',	'CGGA_1402',	'CGGA_1433',	'CGGA_P160',	'CGGA_1262',	'CGGA_1833',	'CGGA_1690',	'CGGA_1815',	'CGGA_P415',	'CGGA_1976',	'CGGA_1498',	'CGGA_1603',	'CGGA_1953',	'CGGA_1462',	'CGGA_1773',	'CGGA_1644',	'CGGA_1255',	'CGGA_1722',	'CGGA_777',	'CGGA_1371',	'CGGA_1663',	'CGGA_1391',	'CGGA_1658',	'CGGA_1682',	'CGGA_1572',	'CGGA_P335',	'CGGA_1478',	'CGGA_1750',	'CGGA_1172',	'CGGA_1807',	'CGGA_1420',	'CGGA_1529',	'CGGA_1542',	'CGGA_1596',	'CGGA_1560',	'CGGA_P112',	'CGGA_P596',	'CGGA_1378',	'CGGA_1728',	'CGGA_P175',	'CGGA_P280',	'CGGA_1410',	'CGGA_1337',	'CGGA_1472',	'CGGA_P136',	'CGGA_1103',	'CGGA_P609',	'CGGA_1985',	'CGGA_1769',	'CGGA_1615',	'CGGA_1564',	'CGGA_1539',	'CGGA_2024',	'CGGA_1735',	'CGGA_P178',	'CGGA_2106',	'CGGA_1041',	'CGGA_2039',	'CGGA_1727',	'CGGA_1605',	'CGGA_2053',	'CGGA_1138',	'CGGA_P99',	'CGGA_139',	'CGGA_1496',	'CGGA_P22',	'CGGA_1425',	'CGGA_1325',	'CGGA_1650',	'CGGA_1164',	'CGGA_1534',	'CGGA_1494',	'CGGA_1955',	'CGGA_1236',	'CGGA_1595',	'CGGA_1870',	'CGGA_P15',	'CGGA_1017',	'CGGA_1611',	'CGGA_1699',	'CGGA_1430',	'CGGA_1571',	'CGGA_1467',	'CGGA_1740',	'CGGA_P154',	'CGGA_1635')
cgga_693_gbm_gsva_wanted <- cgga_693_gbm_gsva[,match(samples_wanted, colnames(cgga_693_gbm_gsva))]
gene_sets_wanted <- c('VERHAAK_GLIOBLASTOMA_MESENCHYMAL', 'PHILIPS_MESENCHYMAL', 'NEFTEL_MES_META_MODULE', 'WANG_mGSC', 'VERHAAK_GLIOBLASTOMA_PRONEURAL', 'NEFTEL_NPC_META_MODULE', 'NEFTEL_OPC_META_MODULE', 'PHILIPS_PRONEURAL', 'WANG_pGSC', 'VERHAAK_GLIOBLASTOMA_CLASSICAL', 'NEFTEL_AC_META_MODULE')
cgga_693_gbm_gsva_wanted <- cgga_693_gbm_gsva_wanted[match(gene_sets_wanted, row.names(cgga_693_gbm_gsva_wanted)),]
my_sample_col <- data.frame(sample = rep(c("CAF_High", "CAF_Low"), c(62,62)))
row.names(my_sample_col) <- colnames(cgga_693_gbm_gsva_wanted)
my_colour = list(sample = c(CAF_High = "red", CAF_Low = "blue"))
pheatmap(cgga_693_gbm_gsva_wanted, scale = "row", cluster_cols = F, cluster_rows = F, color=colorRampPalette(c("blue","white","red"))(50), breaks=seq(-1, 1, length.out = 50),  show_rownames = T, show_colnames = F, annotation_col = my_sample_col, annotation_colors = my_colour)

################################################################################
################################################################################
################################################################################
################################################################################
tcga_gbm_counts <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/TCGA_GBM_counts.txt', header = T, row.names = 1, sep = '\t')
tcga_gbm_counts <- 2^(tcga_gbm_counts)
tcga_gbm_counts <- tcga_gbm_counts -1
tcga_gbm_counts <- round(tcga_gbm_counts)

gene_name_protein_coding <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/PMT_Preliminary_Data/PMT_CAF_accumilation/protein-coding_gene.txt', header = T, row.names = 1, sep = '\t') #protein-coding_gene.txt
tcga_gbm_counts <- tcga_gbm_counts[match(gene_name_protein_coding$symbol, row.names(tcga_gbm_counts)),]
tcga_gbm_counts <- na.omit(tcga_gbm_counts)

tcga_gbm_samples_wanted <- c('TCGA.19.5960.01A',	'TCGA.15.1444.01A',	'TCGA.19.1390.01A',	'TCGA.19.2624.01A',	'TCGA.12.3650.01A',	'TCGA.41.5651.01A',	'TCGA.06.0238.01A',	'TCGA.14.1825.01A',	'TCGA.26.5134.01A',	'TCGA.12.0616.01A',	'TCGA.16.0846.01A',	'TCGA.12.1597.01B',	'TCGA.41.2571.01A',	'TCGA.12.0618.01A',	'TCGA.32.2634.01A',	'TCGA.06.2558.01A',	'TCGA.06.5416.01A',	'TCGA.06.0745.01A',	'TCGA.06.0174.01A',	'TCGA.06.0686.01A',	'TCGA.32.5222.01A',	'TCGA.26.5135.01A',	'TCGA.76.4932.01A',	'TCGA.02.0047.01A',	'TCGA.76.4925.01A',	'TCGA.06.2559.01A',	'TCGA.06.0646.01A',	'TCGA.27.1830.01A',	'TCGA.06.0156.01A',	'TCGA.06.2569.01A',	'TCGA.06.2557.01A',	'TCGA.28.5207.01A',	'TCGA.14.0871.01A',	'TCGA.19.1787.01B',	'TCGA.32.2632.01A',	'TCGA.28.5208.01A',	'TCGA.32.2616.01A',	'TCGA.06.2561.01A',	'TCGA.26.5136.01B',	'TCGA.06.5858.01A',	'TCGA.28.5215.01A',	'TCGA.06.0184.01A',	'TCGA.28.5218.01A',	'TCGA.28.1753.01A',	'TCGA.27.2519.01A',	'TCGA.28.5209.01A',	'TCGA.14.1034.01A',	'TCGA.14.1823.01A',	'TCGA.28.2509.01A',	'TCGA.06.0130.01A',	'TCGA.06.0750.01A',	'TCGA.12.0619.01A',	'TCGA.16.1045.01B',	'TCGA.28.5216.01A',	'TCGA.28.5213.01A',	'TCGA.06.0141.01A',	'TCGA.06.5418.01A',	'TCGA.27.1832.01A',	'TCGA.06.5412.01A',	'TCGA.41.3915.01A',	'TCGA.41.4097.01A',	'TCGA.06.0168.01A',	'TCGA.06.2562.01A',	'TCGA.02.2486.01A',	'TCGA.32.4213.01A',	'TCGA.32.2615.01A',	'TCGA.06.0644.01A',	'TCGA.06.5410.01A',	'TCGA.26.5139.01A',	'TCGA.06.0878.01A',	'TCGA.06.0139.01A',	'TCGA.06.0210.01A',	'TCGA.06.0645.01A',	'TCGA.14.0781.01B',	'TCGA.06.0190.01A',	'TCGA.28.2513.01A',	'TCGA.27.2524.01A',	'TCGA.14.0789.01A',	'TCGA.02.0055.01A')
tcga_gbm_counts_samples_wanted <- tcga_gbm_counts[,match(tcga_gbm_samples_wanted, colnames(tcga_gbm_counts))]
tcga_gbm_counts_samples_wanted <- subset(tcga_gbm_counts_samples_wanted, rowSums(tcga_gbm_counts_samples_wanted)>0)

coldata <- data.frame(row.names = colnames(tcga_gbm_counts_samples_wanted))
coldata$condition <- c(rep("Proneural_CAFs_Low", 15), rep("Proneural_CAFs_High", 14), rep("Mesenchymal_CAFs_Low", 25), rep("Mesenchymal_CAFs_High", 25))

dds <- DESeqDataSetFromMatrix(countData = tcga_gbm_counts_samples_wanted, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

dds2 <- dds[rowSums(counts(dds))>1,]
rld <- vst(dds2,blind=FALSE)
plotPCA <- plotPCA(rld, intgroup = c("condition"))
plotPCA
ggplot(plotPCA$data, aes(plotPCA$data$PC1, plotPCA$data$PC2, col= plotPCA$data$group, fill = plotPCA$data$group))+
  geom_point(shape = 21, col = "black", size = 3) +
  labs(x = "PC1", y = "PC2", title= "TCGA")

plot1 <- ggplot(plotPCA$data, aes(plotPCA$data$PC1, plotPCA$data$PC2, col= plotPCA$data$group, fill = group))+
  geom_point(shape = 21, col = "black", size = 5) +
  labs(x = "PC1", y = "PC2") +
  theme_pubr(legend = c("none")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), title =  element_text(size=14, face="bold", colour = "black"))
  

dens1 <- ggplot(plotPCA$data, aes(plotPCA$data$PC1))+
  geom_density(aes(color = group), size = 1) +
  theme_void() +
  theme(legend.position = "none")

dens1 + plot_spacer() + plot1 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(20, 1), heights = c(1, 1))

normalized_counts <- counts(dds, normalized = T)

rows_wanted <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/PMT_Preliminary_Data/PMT_CAF_accumilation/fulllist.txt', header = T, sep = '\t') #fulllist.txt
normalized_counts_final_wanted <- normalized_counts[match(rows_wanted$Gene, row.names(normalized_counts)),]
normalized_counts_final_wanted <- data.frame(normalized_counts_final_wanted)
dim(normalized_counts_final_wanted)

tcga_gbm_gsva_cell_type <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/GSVA_Analysis/TCGA/TCGA_GBM_gsva.csv', header = T, row.names = 1)
test <- tcga_gbm_gsva_cell_type[match(tcga_gbm_samples_wanted, row.names(tcga_gbm_gsva_cell_type)),]
table(sort(colnames(normalized_counts_final_wanted)==row.names(test)))

my_sample_col <- data.frame(Subtype = factor(rep(c("Proneural_CAF_Low", "Proneural_CAF_High", "Mesenchymal_CAF_Low", "Mesenchymal_CAF_High"), c(15, 14, 25,25))))
row.names(my_sample_col) <- colnames(normalized_counts_final_wanted)
table(sort(row.names(my_sample_col)==row.names(test)))
my_sample_col$CAF <- test$CAF
head(my_sample_col)

my_sample_row <- data.frame(Processes = factor(rep(c("Inflammation", "Invasion", "Mesenchymal", "Proneural"), c(290, 212, 371, 415))))
row.names(my_sample_row) <- row.names(normalized_counts_final_wanted)
head(my_sample_row)

annotation_color <- list(Subtype = c(Proneural_CAF_Low = "#C77CFF" ,Proneural_CAF_High = "#00BFC4" , Mesenchymal_CAF_Low = "#7CAE00" , Mesenchymal_CAF_High = "#F8766D"), CAF = c(colorRampPalette(c("white","grey", "black"))(100)), Processes =c(Invasion = "chocolate", Inflammation = "cadetblue", Mesenchymal = "red", Proneural = "blue"))
pheatmap(normalized_counts_final_wanted, scale = "row", annotation_row = my_sample_row,  annotation_col = my_sample_col, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", color=colorRampPalette(c("blue","white","red"))(50), breaks=seq(-1, 1, length.out = 50), show_rownames = F, show_colnames = F, annotation_colors =  annotation_color, gaps_row=c(290,502,873))

#CGGA_325
CGGA325_gbm_counts <- read.table('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/TCGA_CGGA_DataSets/CGGA.325.gene.count.name.txt', header = T, row.names = 1, sep = '\t')

CGGA325_gbm_counts <- CGGA325_gbm_counts[match(gene_name_protein_coding$symbol, row.names(CGGA325_gbm_counts)),]
CGGA325_gbm_counts <- na.omit(CGGA325_gbm_counts)

CGGA325_gbm_samples_wanted <- c('CGGA_518',	'CGGA_700',	'CGGA_822',	'CGGA_1116',	'CGGA_1060',	'CGGA_747',	'CGGA_1251',	'CGGA_1394',	'CGGA_1375',	'CGGA_1219',	'CGGA_374',	'CGGA_1409',	'CGGA_1450',	'CGGA_J100',	'CGGA_491',	'CGGA_1129',	'CGGA_1170',	'CGGA_632',	'CGGA_899',	'CGGA_D34',	'CGGA_1384',	'CGGA_1197',	'CGGA_1287',	'CGGA_D59',	'CGGA_1095',	'CGGA_719',	'CGGA_761',	'CGGA_499',	'CGGA_658',	'CGGA_1381',	'CGGA_1272',	'CGGA_681',	'CGGA_1475',	'CGGA_710',	'CGGA_1275',	'CGGA_545',	'CGGA_1285',	'CGGA_1136',	'CGGA_679',	'CGGA_1171',	'CGGA_1011',	'CGGA_1053',	'CGGA_495',	'CGGA_309',	'CGGA_1234',	'CGGA_782',	'CGGA_1074',	'CGGA_272',	'CGGA_1218',	'CGGA_483',	'CGGA_789',	'CGGA_D02',	'CGGA_1019',	'CGGA_902',	'CGGA_1001',	'CGGA_1114',	'CGGA_1105',	'CGGA_1224',	'CGGA_1083',	'CGGA_525',	'CGGA_1015',	'CGGA_876',	'CGGA_D37',	'CGGA_1180',	'CGGA_D51',	'CGGA_604',	'CGGA_1313',	'CGGA_1077',	'CGGA_1049',	'CGGA_D35',	'CGGA_842',	'CGGA_D36',	'CGGA_1338',	'CGGA_D26',	'CGGA_624',	'CGGA_802',	'CGGA_1124',	'CGGA_1240',	'CGGA_1342',	'CGGA_1045',	'CGGA_1177',	'CGGA_859',	'CGGA_1035',	'CGGA_318')
CGGA325_gbm_counts_samples_wanted <- CGGA325_gbm_counts[,match(CGGA325_gbm_samples_wanted, colnames(CGGA325_gbm_counts))]
CGGA325_gbm_counts_samples_wanted <- subset(CGGA325_gbm_counts_samples_wanted, rowSums(CGGA325_gbm_counts_samples_wanted)>0)

coldata <- data.frame(row.names = colnames(CGGA325_gbm_counts_samples_wanted))
coldata$condition <- c(rep("Proneural_CAFs_Low", 17), rep("Proneural_CAFs_High", 16), rep("Mesenchymal_CAFs_Low", 26), rep("Mesenchymal_CAFs_High", 25))

dds <- DESeqDataSetFromMatrix(countData = CGGA325_gbm_counts_samples_wanted, colData = coldata, design = ~ condition)
dds <- DESeq(dds)

dds2 <- dds[rowSums(counts(dds))>2,]
rld <- vst(dds2,blind=FALSE)
plotPCA <- plotPCA(rld, intgroup = c("condition"))
plotPCA
ggplot(plotPCA$data, aes(plotPCA$data$PC1, plotPCA$data$PC2, col= plotPCA$data$group, fill = plotPCA$data$group))+
  geom_point(shape = 21, col = "black", size = 5) +
  labs(x = "PC1", y = "PC2", title= "CGGA_325")

plot1 <- ggplot(plotPCA$data, aes(plotPCA$data$PC1, plotPCA$data$PC2, col= plotPCA$data$group, fill = group))+
  geom_point(shape = 21, col = "black", size = 5) +
  labs(x = "PC1", y = "PC2") +
  theme_pubr(legend = c("none")) +
  theme(axis.text.x = element_text(angle = 0, vjust = 1, hjust=1, face = "bold", size = 12), axis.text.y = element_text(face = "bold", size = 12), axis.title.x = element_text(size=14, face="bold", colour = "black"), axis.title.y = element_text(size=14, face="bold", colour = "black"), title =  element_text(size=14, face="bold", colour = "black"))


dens1 <- ggplot(plotPCA$data, aes(plotPCA$data$PC1))+
  geom_density(aes(color = group), size = 1) +
  theme_void() +
  theme(legend.position = "none")

dens1 + plot_spacer() + plot1 + 
  plot_layout(ncol = 2, nrow = 2, widths = c(20, 1), heights = c(1, 1))

normalized_counts <- counts(dds, normalized = T)

normalized_counts_final_wanted <- normalized_counts[match(rows_wanted$Gene, row.names(normalized_counts)),]
normalized_counts_final_wanted <- data.frame(normalized_counts_final_wanted)
dim(normalized_counts_final_wanted)

CGGA_325_GBM_GSVA_celltype <- read.csv('/Users/phillipgalbo/Desktop/GBM_CAF_Paper/GSVA_Analysis/CGGA_325/cgga_325_gsva.csv', header = T, row.names = 1)
test <- CGGA_325_GBM_GSVA_celltype[match(CGGA325_gbm_samples_wanted, row.names(CGGA_325_GBM_GSVA_celltype)),]
table(sort(colnames(normalized_counts_final_wanted)==row.names(test)))

my_sample_col <- data.frame(Subtype = factor(rep(c("Proneural_CAF_Low", "Proneural_CAF_High", "Mesenchymal_CAF_Low", "Mesenchymal_CAF_High"), c(17, 16, 26,25))))
row.names(my_sample_col) <- colnames(normalized_counts_final_wanted)
table(sort(row.names(my_sample_col)==row.names(test)))
my_sample_col$CAF <- test$CAF
head(my_sample_col)

my_sample_row <- data.frame(Processes = factor(rep(c("Inflammation", "Invasion", "Mesenchymal", "Proneural"), c(290, 212, 371, 415))))
row.names(my_sample_row) <- row.names(normalized_counts_final_wanted)
head(my_sample_row)

annotation_color <- list(Subtype = c(Proneural_CAF_Low = "#C77CFF" ,Proneural_CAF_High = "#00BFC4" , Mesenchymal_CAF_Low = "#7CAE00" , Mesenchymal_CAF_High = "#F8766D"), CAF = c(colorRampPalette(c("white","grey", "black"))(100)), Processes =c(Invasion = "chocolate", Inflammation = "cadetblue", Mesenchymal = "burlywood", Proneural = "darkseagreen"))
pheatmap(normalized_counts_final_wanted, scale = "row", annotation_row = my_sample_row,  annotation_col = my_sample_col, cluster_rows = F, cluster_cols = T, clustering_method = "ward.D2", color=colorRampPalette(c("blue","white","red"))(50), breaks=seq(-1, 1, length.out = 50), show_rownames = F, show_colnames = F, annotation_colors =  annotation_color, gaps_row=c(290,502,873))























