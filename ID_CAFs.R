library(devtools)
library(RaceID)
library(biomaRt)
library(scran)
library(ggplot2)
library(ggrepel)

#Load the required files
gbm_seq <- read.table('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/normal_cells_with_cafs.txt', header = T, row.names = 1, sep = "\t")
sc <- readRDS('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/sc_rfcorrect.rds')
sc_umap <- readRDS('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/sc_umap.rds')
umap_coordiantes <- read.table(file.choose(), header = T, row.names = 1, sep = '\t') #read: umap_coordinates.txt which is located on Deyous computer via the following path: /home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235
ggplot(umap_coordiantes, aes(V1, V2)) + geom_point() 

vascular_cells <- subset(umap_coordiantes, umap_coordiantes$V1 > -5 & umap_coordiantes$V1 < 5 & umap_coordiantes$V2 > -5 & umap_coordiantes$V2 < -2.5)
ggplot(vascular_cells, aes(V1, V2)) + geom_point() + xlim(-10, 18) + ylim(-12, 12)

cluster_want_one <- subset(umap_coordiantes, umap_coordiantes$V1 > -3 & umap_coordiantes$V1 < 1 & umap_coordiantes$V2 > 0 & umap_coordiantes$V2 < 5)
ggplot(cluster_want_one, aes(V1, V2)) + geom_point() + xlim(-10, 18) + ylim(-12, 12)

cluster_want_two <- subset(umap_coordiantes, umap_coordiantes$V1 > -5 & umap_coordiantes$V1 < -2 & umap_coordiantes$V2 > 2 & umap_coordiantes$V2 < 5)
ggplot(cluster_want_two, aes(V1, V2)) + geom_point() + xlim(-10, 18) + ylim(-12, 12)


#Code to run on Deyou's Cluster
gbm_seq <- read.table('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/normal_cells_with_cafs.txt', header = T, row.names = 1, sep = "\t")
sc <- readRDS('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/sc_rfcorrect.rds')
sc_umap <- readRDS('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/sc_umap.rds')

#Make UMAP with CAFs
#Get Random Clusters
vascular_cells <- subset(sc_umap@umap, sc_umap@umap$V1 > -5 & sc_umap@umap$V1 < 5 & sc_umap@umap$V2 > -5 & sc_umap@umap$V2 < -2.5)
cluster_want_one <- subset(sc_umap@umap, sc_umap@umap$V1 > -3 & sc_umap@umap$V1 < 1 & sc_umap@umap$V2 > 0 & sc_umap@umap$V2 < 5)
cluster_want_two <- subset(sc_umap@umap, sc_umap@umap$V1 > -5 & sc_umap@umap$V1 < -2 & sc_umap@umap$V2 > 2 & sc_umap@umap$V2 < 5)

#Get CAFs
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

#Add other cells to cparts as own cluster
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
pdf("/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/ggplot_of_clusters_with_various_clusters_labeled.pdf")
ggplot(sc_umap@umap, aes(x=sc_umap@umap$V1, y=sc_umap@umap$V2, color = test$test)) + geom_point() + labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") + theme(legend.title.align = 0.5) 
dev.off()

cdiff_172 <- clustdiffgenes(sc, cl = 172, pvalue = 0.0499)

#Make UMAP without CAFs
vascular_cells <- subset(sc_umap@umap, sc_umap@umap$V1 > -5 & sc_umap@umap$V1 < 5 & sc_umap@umap$V2 > -5 & sc_umap@umap$V2 < -2.5)
cluster_want_one <- subset(sc_umap@umap, sc_umap@umap$V1 > -3 & sc_umap@umap$V1 < 1 & sc_umap@umap$V2 > 0 & sc_umap@umap$V2 < 5)
cluster_want_two <- subset(sc_umap@umap, sc_umap@umap$V1 > -5 & sc_umap@umap$V1 < -2 & sc_umap@umap$V2 > 2 & sc_umap@umap$V2 < 5)

interest <- sample(row.names(cluster_want_one), 119)
for(i in interest){
  sc_umap@cluster$kpart[names(sc_umap@cluster$kpart) == i] = 174
}

interest <- sample(row.names(cluster_want_two), 60)
for(i in interest){
  sc_umap@cluster$kpart[names(sc_umap@cluster$kpart) == i] = 175
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
test$test <- factor(test$test, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "Cluster 11"))
pdf("/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/ggplot_of_clusters_with_various_clusters_labeled_without_cafs.pdf")
ggplot(sc_umap@umap, aes(x=sc_umap@umap$V1, y=sc_umap@umap$V2, color = test$test)) + geom_point() + labs(x = "UMAP 1", y = "UMAP 2", colour = "Cluster") + theme(legend.title.align = 0.5) 
dev.off()

#Generate CIBERSORTx matrix
test <- sc_umap@cluster$kpart
test <- data.frame(test)

cluster1_names <- subset(test, test == 1)
gbm_seq_cluster1 <- gbm_seq [,match(row.names(cluster1_names), colnames(gbm_seq))]
dim(gbm_seq_cluster1)
gbm_seq_cluster1[1:5,1:5]
col.data.1 <- c("Cell Cluster 1")
colnames(gbm_seq_cluster1)[1:1912] <- col.data.1
gbm_seq_cluster1[1:5,1:5]

cluster2_names <- subset(test, test == 2)
gbm_seq_cluster2 <- gbm_seq [,match(row.names(cluster2_names), colnames(gbm_seq))]
dim(gbm_seq_cluster2)
gbm_seq_cluster2[1:5,1:5]
col.data.2 <- c("Cell Cluster 2")
colnames(gbm_seq_cluster2)[1:837] <- col.data.2
gbm_seq_cluster2[1:5,1:5]

cluster3_names <- subset(test, test == 3)
gbm_seq_cluster3 <- gbm_seq [,match(row.names(cluster3_names), colnames(gbm_seq))]
dim(gbm_seq_cluster3)
gbm_seq_cluster3[1:5,1:5]
col.data.3 <- c("Cell Cluster 3")
colnames(gbm_seq_cluster3)[1:1197] <- col.data.3
gbm_seq_cluster3[1:5,1:5]

cluster4_names <- subset(test, test == 4)
gbm_seq_cluster4 <- gbm_seq [,match(row.names(cluster4_names), colnames(gbm_seq))]
dim(gbm_seq_cluster4)
gbm_seq_cluster4[1:5,1:5]
col.data.4 <- c("Cell Cluster 4")
colnames(gbm_seq_cluster4)[1:281] <- col.data.4
gbm_seq_cluster4[1:5,1:5]

cluster5_names <- subset(test, test == 5)
gbm_seq_cluster5 <- gbm_seq [,match(row.names(cluster5_names), colnames(gbm_seq))]
dim(gbm_seq_cluster5)
gbm_seq_cluster5[1:5,1:5]
col.data.5 <- c("Cell Cluster 5")
colnames(gbm_seq_cluster5)[1:1986] <- col.data.5
gbm_seq_cluster5[1:5,1:5]

cluster6_names <- subset(test, test == 6)
gbm_seq_cluster6 <- gbm_seq [,match(row.names(cluster6_names), colnames(gbm_seq))]
dim(gbm_seq_cluster6)
gbm_seq_cluster6[1:5,1:5]
col.data.6 <- c("Cell Cluster 6")
colnames(gbm_seq_cluster6)[1:1324] <- col.data.6
gbm_seq_cluster6[1:5,1:5]

cluster7_names <- subset(test, test == 7)
gbm_seq_cluster7 <- gbm_seq [,match(row.names(cluster7_names), colnames(gbm_seq))]
dim(gbm_seq_cluster7)
gbm_seq_cluster7[1:5,1:5]
col.data.7 <- c("Cell Cluster 7")
colnames(gbm_seq_cluster7)[1:1867] <- col.data.7
gbm_seq_cluster7[1:5,1:5]

cluster8_names <- subset(test, test == 8)
gbm_seq_cluster8 <- gbm_seq [,match(row.names(cluster8_names), colnames(gbm_seq))]
dim(gbm_seq_cluster8)
gbm_seq_cluster8[1:5,1:5]
col.data.8 <- c("Cell Cluster 8")
colnames(gbm_seq_cluster8)[1:48] <- col.data.8
gbm_seq_cluster8[1:5,1:5]

cluster9_names <- subset(test, test == 172) #CAFs
gbm_seq_cluster9 <- gbm_seq [,match(row.names(cluster9_names), colnames(gbm_seq))]
dim(gbm_seq_cluster9)
gbm_seq_cluster9[1:5,1:5]
col.data.9 <- c("Cell Cluster 9")
colnames(gbm_seq_cluster9)[1:265] <- col.data.9
gbm_seq_cluster9[1:5,1:5]

cluster10_names <- subset(test, test == 173)
gbm_seq_cluster10 <- gbm_seq [,match(row.names(cluster10_names), colnames(gbm_seq))]
dim(gbm_seq_cluster10)
gbm_seq_cluster10[1:5,1:5]
col.data.10 <- c("Cell Cluster 10")
colnames(gbm_seq_cluster10)[1:115] <- col.data.10
gbm_seq_cluster10[1:5,1:5]

cluster11_names <- subset(test, test == 174)
gbm_seq_cluster11 <- gbm_seq [,match(row.names(cluster11_names), colnames(gbm_seq))]
dim(gbm_seq_cluster11)
gbm_seq_cluster11[1:5,1:5]
col.data.11 <- c("Cell Cluster 11")
colnames(gbm_seq_cluster11)[1:117] <- col.data.11
gbm_seq_cluster11[1:5,1:5]

cluster12_names <- subset(test, test == 175)
gbm_seq_cluster12 <- gbm_seq [,match(row.names(cluster12_names), colnames(gbm_seq))]
dim(gbm_seq_cluster12)
gbm_seq_cluster12[1:5,1:5]
col.data.12 <- c("Cell Cluster 12")
colnames(gbm_seq_cluster12)[1:51] <- col.data.12
gbm_seq_cluster12[1:5,1:5]

table(row.names(gbm_seq_cluster1)==row.names(gbm_seq_cluster2))
table(row.names(gbm_seq_cluster3)==row.names(gbm_seq_cluster2))
table(row.names(gbm_seq_cluster3)==row.names(gbm_seq_cluster4))
table(row.names(gbm_seq_cluster4)==row.names(gbm_seq_cluster5))
table(row.names(gbm_seq_cluster5)==row.names(gbm_seq_cluster6))
table(row.names(gbm_seq_cluster6)==row.names(gbm_seq_cluster7))
table(row.names(gbm_seq_cluster8)==row.names(gbm_seq_cluster9))
table(row.names(gbm_seq_cluster9)==row.names(gbm_seq_cluster10))
table(row.names(gbm_seq_cluster10)==row.names(gbm_seq_cluster11))
table(row.names(gbm_seq_cluster11)==row.names(gbm_seq_cluster12))

gbm_seq_race_clusters <- cbind(gbm_seq_cluster1, gbm_seq_cluster2, gbm_seq_cluster3, gbm_seq_cluster4, gbm_seq_cluster5, gbm_seq_cluster6, gbm_seq_cluster7, gbm_seq_cluster8, gbm_seq_cluster9, gbm_seq_cluster10, gbm_seq_cluster11, gbm_seq_cluster12)
write.table(gbm_seq_race_clusters, '/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/cibersortx_file_with_random_clusters_labeled.txt', col.names = T, row.names = T, sep = "\t", quote = F)
cibersortx_file <- read.table('/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/cibersortx_file_with_random_clusters_labeled.txt', header = T, row.names = 1, sep = "\t")

write.table(test, '/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/umap_cluster_kpart.txt', col.names = T, row.names = T, sep = "\t", quote = F)

#Get marker genes for CAF cluster
interest <- sample(row.names(cluster_want_one), 119)
for(i in interest){
  sc@cpart[names(sc@cpart) == i] = 174
}

interest <- sample(row.names(cluster_want_two), 60)
for(i in interest){
  sc@cpart[names(sc@cpart) == i] = 175
}

interest <- sample(row.names(caf_like_cells_t_wanted), 278)
for(i in interest){
  sc@cpart[names(sc@cpart) == i] = 172
}

interest <- sample(row.names(vascular_cells), 115)
for(i in interest){
  sc@cpart[names(sc@cpart) == i] = 173
}

cdiff_172 <- clustdiffgenes(sc, cl = 172, pvalue = 0.0499)
write.table(cdiff_172, file = '/home/phillip/GBM_sc_RACE_23K_analysis/race_analysis_half_cells/first_half/cafs_235/CAFcluster_MarkerGenes.csv', col.names = T, row.names = T, quote = F)







