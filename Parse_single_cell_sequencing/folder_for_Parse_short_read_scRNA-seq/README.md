## Making Markdown for scRNA-seq

The markdown files and scripts are currnetly located in the Biowulf data folder `parse_single_cell/R_scripts`

- parse install
- how to run
- the output
- Seurat
- cell annotation
- public database


-------
### if wnat, can re-do clustering and homony to see test some para to see good clustering, copy from Jyieon's grop clustering

```R

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
library(SeuratData)

lung_4 <- readRDS(file = "lung_merge_all.rds")

dim(lung_4)

table(lung_4@meta.data$orig.ident)

####Harmony on RNA data

DefaultAssay(lung_4) <- "RNA"

lung_4 <- SCTransform(lung_4, verbose = TRUE)

##Set default assay to SCT

DefaultAssay(lung_4) <- "SCT"

#dimensional reduction and clustering

lung_4 <- ScaleData(lung_4, verbose = TRUE)
lung_4 <- RunPCA(lung_4, npcs = 50, verbose = TRUE)
lung_4 <- RunUMAP(lung_4, reduction = "pca", dims = 1:50, reduction.name = 'umap.RNA', reduction.key = 'rnaUMAP_')
lung_4 <- FindNeighbors(lung_4, reduction = "pca", dims = 1:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p1 <- DimPlot(lung_4, reduction = "umap.RNA", group.by = "orig.ident")
p2 <- DimPlot(lung_4, reduction = "umap.RNA", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_SCT.pdf", p1+p2, width = 24, height = 12)

#correcting batch effect

lung_4 <- RunHarmony(object = lung_4, group.by.vars = 'orig.ident',
reduction = 'pca',
assay.use = 'SCT',
reduction.save = "harmony_SCT",
project.dim = FALSE)

lung_4 <- RunUMAP(lung_4, dims = 1:50, reduction = 'harmony_SCT', reduction.name = "umap.SCT_harm", reduction.key = "sctUMAPharm_")
lung_4 <- FindNeighbors(lung_4, reduction = "harmony_SCT", dims = 1:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p3 <- DimPlot(lung_4, reduction = "umap.SCT_harm", group.by = "orig.ident")
p4 <- DimPlot(lung_4, reduction = "umap.SCT_harm", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_SCT_Har.pdf", p3+p4, width = 24, height = 12)

####Harmony on ATAC data

#dimensional reduction and clustering

DefaultAssay(lung_4) <- "ATAC"

lung_4 <- RunTFIDF(lung_4)
lung_4 <- FindTopFeatures(lung_4, min.cutoff = 'q0')
lung_4 <- RunSVD(lung_4)
lung_4 <- RunUMAP(lung_4, reduction = 'lsi', dims = 2:50, reduction.name = "umap.ATAC", reduction.key = "atacUMAP_")

lung_4 <- FindNeighbors(lung_4, reduction = "lsi", dims = 2:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p5 <- DimPlot(lung_4, reduction = "umap.ATAC", group.by = "orig.ident")
p6 <- DimPlot(lung_4, reduction = "umap.ATAC", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_ATAC.pdf", p5+p6, width = 24, height = 12)

#correcting batch effect

lung_4 <- RunHarmony(object = lung_4, group.by.vars = 'orig.ident',
reduction = 'lsi',
assay.use = 'ATAC',
reduction.save = "harmony_ATAC",
project.dim = FALSE)

lung_4 <- RunUMAP(lung_4, dims = 2:50, reduction = 'harmony_ATAC', reduction.name = "umap.ATAC_harm", reduction.key = "atacUMAPharm_")
lung_4 <- FindNeighbors(lung_4, reduction = "harmony_ATAC", dims = 2:50)
lung_4 <- FindClusters(lung_4, resolution = 0.5)

p7 <- DimPlot(lung_4, reduction = "umap.ATAC_harm", group.by = "orig.ident")
p8 <- DimPlot(lung_4, reduction = "umap.ATAC_harm", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_ATAC_Har.pdf", p7+p8, width = 24, height = 12)


####WNN on two modalities after harmony

lung_4 <- FindMultiModalNeighbors(lung_4, reduction.list = list("harmony_SCT", "harmony_ATAC"), dims.list = list(1:50, 2:50))
lung_4 <- RunUMAP(lung_4, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")

#final clustering

lung_4 <- FindClusters(lung_4, graph.name = "wsnn", algorithm = 3, resolution = 0.5)

p9 <- DimPlot(lung_4, reduction = "wnn.umap", group.by = "orig.ident")
p10 <- DimPlot(lung_4, reduction = "wnn.umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE)

ggsave("afterintegrate_UMAP_WNN_Har.pdf", p9+p10, width = 24, height = 12)

saveRDS(lung_4, file = "lung_integrate_all.rds")





```




### Additioanl analysis

**The Code from Jiyeon's lab for sub-clustering is good for me to sub-cluster our Parse-SC bladder tissue bladder urithrial group. **

```
####################################################################################
# Date: 04/01/2025
# Title: Sub-cluster
# Info: The script subset the target cluster, and re-do clustering etc to get
# the more detealied info cluster in that cluster
#####################################################################################

library(Seurat)
library(Signac)
library(harmony)
library(ggplot2)


setwd('/data/Choi_lung/ChiaHan/')



lung_4 <- readRDS("CHL_pbmc_integrate_all_removal_combined_res1_annot_with_peak_chranno_TF_FOOTPRINTED_Linkage1MB.rds")
# set cell 
cell_of_interest <- c("NK",'T','Macrophage','Monocyte','Lymphatic','Artery','Vein','Capillary','AT1','AT2','Club',
                      'Ciliated','Goblet','Basal','Fibroblast','SMC')

cell_of_interest_BLUE_Immune <- c("NK",'T','Macrophage','Monocyte')
cell_of_interest_Red_endothelial <- c('Lymphatic','Artery','Vein','Capillary')
cell_of_interest_Green_epithelial <- c('AT1','AT2','Club','Ciliated')  
cell_of_interest_Yellow_stromal <- c('Goblet','Basal','Fibroblast','SMC')  

lung_4$ite <- as.character(lung_4$CellType)

Idents(lung_4) <- "CellType"

for (i in cell_of_interest){
  

  
# the original code:  lambda = 0.5,  Default theta = 2, sigma = 0.3

#### The new set #####  
# NK: 0.5,2,0.4
# For T cell : lambda = 0.8,Default theta = 2,sigma = 0.3
# Monocyte: 0.8,2,0.5
# Artery: 0.8,2,0.3
# Vien: 0.8,2,0.2
# Club: 0.8,2,0.3
# Basal: 0.1,2,0.9 
# Goblet: 0.1,2,0.6
# Fibroblast: 0.1,2,0.3; also all the res in FindMarker is 0.015
# lymphatic: 0.4,2,0.4

input_lambda <- 0.4
input_theta <- 2
input_sigma <- 0.4



#i <- cell_of_interest[5]

pbmc_0 <- subset(lung_4, idents = i)

DefaultAssay(pbmc_0) <- "SCT"

pbmc_0 <- ScaleData(pbmc_0, verbose = TRUE)
pbmc_0 <- RunPCA(pbmc_0, npcs = 50, verbose = TRUE) #进行pca降维，维度选取为50
#非线性降维umap，对pca降维的结果进行继续降维，维度为50，结果为uamp.RNA存在reduction中
pbmc_0 <- RunUMAP(pbmc_0, reduction = "pca", dims = 1:50, reduction.name = 'umap.RNA')
#聚类
pbmc_0 <- FindNeighbors(pbmc_0, reduction = "pca", dims = 1:50)
#  the res to change or not?
pbmc_0 <- FindClusters(pbmc_0, resolution = 0.1)



# pbmc_0 <- RunHarmony(object = pbmc_0, group.by.vars = 'orig.ident', lambda = 0.5, sigma = 0.3,
#                      reduction = 'pca',
#                      assay.use = 'SCT',
#                      reduction.save = "harmony_SCT",
#                      project.dim = FALSE)

pbmc_0 <- RunHarmony(object = pbmc_0, group.by.vars = 'orig.ident', lambda = input_lambda, sigma = input_sigma,
                    
                     theta=input_theta,
                     reduction = 'pca',
                     assay.use = 'SCT',
                     reduction.save = "harmony_SCT",
                     project.dim = FALSE)

#对harmony之后的pca,进行umap降维，存为umap.SCT_harm,reduction.key是列名
pbmc_0 <- RunUMAP(pbmc_0, dims = 1:50, reduction = 'harmony_SCT', reduction.name = "umap.SCT_harm", reduction.key = "sctUMAPharm_")
pbmc_0 <- FindNeighbors(pbmc_0, reduction = "harmony_SCT", dims = 1:50)
pbmc_0 <- FindClusters(pbmc_0, resolution = 0.1)


DefaultAssay(pbmc_0) <- "ATAC"

pbmc_0 <- RunTFIDF(pbmc_0) #normalization
#类似于RNA数据里的选择变化较大的数据,找到出现频率高的feature
pbmc_0 <- FindTopFeatures(pbmc_0, min.cutoff = 'q0')
#降维
pbmc_0 <- RunSVD(pbmc_0)
pbmc_0 <- RunUMAP(pbmc_0, reduction = 'lsi', dims = 2:50, reduction.name = "umap.ATAC", reduction.key = "atacUMAP_")

pbmc_0 <- FindNeighbors(pbmc_0, reduction = "lsi", dims = 2:50)

pbmc_0 <- FindClusters(pbmc_0, resolution = 0.5)





# pbmc_0 <- RunHarmony(object = pbmc_0, group.by.vars = 'orig.ident', lambda = 0.5, sigma = 0.3,
#                      reduction = 'lsi',
#                      assay.use = 'ATAC',
#                      reduction.save = "harmony_ATAC",
#                      project.dim = FALSE)


pbmc_0 <- RunHarmony(object = pbmc_0, group.by.vars = 'orig.ident', lambda = input_lambda, sigma = input_sigma,
                     theta= input_theta,
                     reduction = 'lsi',
                     assay.use = 'ATAC',
                     reduction.save = "harmony_ATAC",
                     project.dim = FALSE)

pbmc_0 <- RunUMAP(pbmc_0, dims = 2:50, reduction = 'harmony_ATAC', reduction.name = "umap.ATAC_harm", reduction.key = "atacUMAPharm_")
pbmc_0 <- FindNeighbors(pbmc_0, reduction = "harmony_ATAC", dims = 2:50)

pbmc_0 <- FindClusters(pbmc_0, resolution = 0.1)


pbmc_0 <- FindMultiModalNeighbors(pbmc_0, reduction.list = list("harmony_SCT", "harmony_ATAC"), dims.list = list(1:50, 2:50))
pbmc_0 <- RunUMAP(pbmc_0, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
pbmc_0 <- FindClusters(pbmc_0, graph.name = "wsnn", algorithm = 3, resolution = 0.05)

pbmc_0$ite <- as.character(pbmc_0$seurat_clusters)

for (j in levels(pbmc_0$seurat_clusters)) {
  pbmc_0$ite[pbmc_0$ite == j] <- paste0(i,"_", j)
}

# p9 <- DimPlot(pbmc_0, reduction = "wnn.umap", group.by = "orig.ident", raster = FALSE)

p10  <- DimPlot(pbmc_0, reduction = "wnn.umap", group.by = "ite", 
               label = FALSE, repel = TRUE, raster = FALSE) +ggtitle(paste0(i))

#ggsave(paste0('Fig_S18_Sub_cluster_figure_updated/',i,'_wnn_0.8_2_0.3.pdf'), p10, width = 5.6, height = 5)

print('save!')

rm(pbmc_0)
gc()

}
 

### Make datafreame

#df_all <- data.frame()

# Get the frequency table
tbl <- table(pbmc_0@meta.data$ite)

# Convert to a data frame with the desired column names
df <- data.frame(
  subcluster = names(tbl),
  number     = as.numeric(tbl)
)

df_all <- rbind(df,df_all)

# DotPlot 
DefaultAssay(pbmc_0) <- 'RNA'

feature_gene_T <- c('CD4','CD40LG','TNFRSF25','CD28','TRAT1','CTLA4','CD8A','CD8B','TRGC2')
feature_gene_Fibroblast <- c('MFAP5' ,'SCARA5', 'PI16' ,'SPINT2', 'LIMCH1', 'FGFR4')

pgene <- DotPlot(pbmc_0,features = feature_gene_Fibroblast,group.by = "ite") + xlab(NULL) +  ylab(NULL) + theme(axis.text.x = element_text(angle = 45, hjust = 1))  
# ggsave(paste0('Fig_S18_Sub_cluster_figure/',i,'_FigS19_marker.pdf'), pgene, width = 8, height = 6)

```

#### Running cellchat for detecting cell-cell communt strength?

Cellchat Github: 

https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/Comparison_analysis_of_multiple_datasets.html

```{R}

### Running Cellchat

library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(dplyr)
library(ggplot2)
library(sctransform)
library(patchwork)
library(GenomeInfoDb)
library(tidyverse)
library(grid)
library(readr)
library(harmony)
library(GenomicRanges)
# library(SeuratData)


# load RDS 
setwd('/data/Choi_lung/ChiaHan/')

##############################################################################
# Description: CellCHAT
##############################################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(SingleR)
library(celldex)
library(cowplot)
library(tidyverse)
library(CellChat)


lung_ns <- readRDS('/data/leec20/Jiyeon_lab_single_cell_quest/Scripts/Tmp_RDS/lung_ns.rds') 
lung_s <- readRDS('/data/leec20/Jiyeon_lab_single_cell_quest/Scripts/Tmp_RDS/lung_S.rds') 


#create cellchat object

DefaultAssay(lung_ns) <- "RNA"
DefaultAssay(lung_s) <- "RNA"

Idents(lung_ns) <- "CellType"
Idents(lung_s) <- "CellType"


# Normalized RNA 
lung_ns <- NormalizeData(lung_ns)
lung_s <-  NormalizeData(lung_s)

print('done nomralized')

labels_Non <- Idents(lung_ns)
labels_SM <- Idents(lung_s)

lung_ns@meta.data$samples <- lung_ns@meta.data$SampleID
cellchat_ns <- createCellChat(lung_ns, meta = lung_ns@meta.data, group.by = "CellType",assay = "RNA")

lung_s@meta.data$samples <- lung_s@meta.data$SampleID
cellchat_s <- createCellChat(lung_s, meta = lung_s@meta.data, group.by = "CellType",assay = "RNA")

cellchat_ns@DB <- CellChatDB.human
cellchat_s@DB <- CellChatDB.human

#calculate the communication probability

cellchat_ns <- subsetData(cellchat_ns)
cellchat_ns <- identifyOverExpressedGenes(cellchat_ns)
cellchat_ns <- identifyOverExpressedInteractions(cellchat_ns)
cellchat_ns <- computeCommunProb(cellchat_ns)
cellchat_ns <- computeCommunProbPathway(cellchat_ns)
cellchat_ns <- aggregateNet(cellchat_ns)

cellchat_s <- subsetData(cellchat_s)
cellchat_s <- identifyOverExpressedGenes(cellchat_s)
cellchat_s <- identifyOverExpressedInteractions(cellchat_s)
cellchat_s <- computeCommunProb(cellchat_s)
cellchat_s <- computeCommunProbPathway(cellchat_s)
cellchat_s <- aggregateNet(cellchat_s)

cellchat_ns <- netAnalysis_computeCentrality(object = cellchat_ns, slot.name = "netP")
cellchat_s <- netAnalysis_computeCentrality(object = cellchat_s, slot.name = "netP")

object.list <- list(ns = cellchat_ns, s = cellchat_s)

cellchat_com <- mergeCellChat(object.list, add.names = names(object.list))

print('Celled')

# save rds

saveRDS(cellchat_com,'Cellchat_result/cellchat_result_list.rds')
print('save!')

#save the probability
ns.netp <- subsetCommunication(cellchat_ns, slot.name = 'netP')
write.csv(ns.netp, "Cellchat_result/pathway_nonsmoke_CHL_v2.csv")

ns.net <- subsetCommunication(cellchat_ns, slot.name = 'net')
write.csv(ns.net, "Cellchat_result/L_R_nonsmoke_CHL_v2.csv")

### DONE ####
print('save!')


### Plotting Cellchat
# CellChat plotting
# Date: 07142025
# CHL, and thanks Bolun for helping a lot
# load files 
library(ggplot2)
library(dplyr)
library(CellChat)
library(tidyr)
library(pheatmap)
library(reshape2)
setwd("/data/Choi_lung/ChiaHan/")
# laod dataset 
cellchat_result_list <- readRDS('Cellchat_result/cellchat_result_list.rds')
pathway_smoke_CHL_v2 <- read.csv('Cellchat_result/pathway_smoke_CHL_v2.csv')
pathway_nonsmoke_CHL_v2 <- read.csv('Cellchat_result/pathway_nonsmoke_CHL_v2.csv')

##### Summary table

pathway_nonsmoke_CHL_v2$logP <- log(pathway_nonsmoke_CHL_v2$prob)
pathway_smoke_CHL_v2$logP <- log(pathway_smoke_CHL_v2$prob)

tmp_ns_v2 <- pathway_nonsmoke_CHL_v2 %>% group_by(pathway_name) %>% summarise(avg = exp(mean(logP)),
                                                                              mean = mean(prob),
                                                                              sum = sum(prob))

tmp_s_v2 <- pathway_smoke_CHL_v2 %>% group_by(pathway_name) %>% summarise(avg = exp(mean(logP)),
                                                                          mean = mean(prob),
                                                                          sum = sum(prob))

tmp <- left_join(tmp_ns_v2 , tmp_s_v2, by = "pathway_name", suffix = c(".ns", ".s"))

tmp$mean.s[is.na(tmp$mean.s)] <- 0
tmp$avg.s[is.na(tmp$avg.s)] <- 0
tmp$sum.s[is.na(tmp$sum.s)] <- 0
sum(is.na(tmp))
tmp$delta_mean <- tmp$mean.s - tmp$mean.ns
tmp$delta_avg <- tmp$avg.s - tmp$avg.ns
tmp$delta_sum <- tmp$sum.s - tmp$sum.ns
# Slect top 5 for sk and non-sk based on delta_sum wiht p-val is sig for heatmap 

rankdf <- rankNet(cellchat_result_list, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = T, do.stat = TRUE,return.data = T)
rankdf <- rankdf$signaling.contribution
rankdf <- rankdf %>% filter(pvalues < 0.05)
rankdf$relative_sum <- rankdf$contribution.scaled/max(rankdf$contribution.scaled)
pathway_list <- unique(rankdf$name)

# since compare rankdf the "LAMININ" is not sig, so chose top 6 and remove LAMININ
tmp <- subset(tmp, pathway_name %in% pathway_list)
top5 <- head(tmp[order(tmp$delta_sum),"pathway_name"], 5)
bottom5 <- head(tmp[order(tmp$delta_sum, decreasing = TRUE),"pathway_name"], 5)


# so gene of interest 

GeneInterest <- c(top5$pathway_name,bottom5$pathway_name)
#write.table(tmp, file = "/data/Choi_lung/ChiaHan/Cellchat_result/CellChat_sum.txt", sep = "\t", row.names = FALSE, quote = FALSE)


############################
#####  NOW we plot #########
############################

# chage the order 
cell_order <- c('Lymphatic','Artery','Vein',	'Capillary'	,'Fibroblast','SMC',
                'Mesothelial',	'MyoFib',	'AT2'	,'AT1'	,'Club',	'Ciliated',
                'Goblet',	'Basal',	'AT1_AT2',	'AT2_pro',	'NK',	'T',	'Macrophage',
                'Monocyte',	'B',	'DC',	'NK_T')


################################
# plot set1 and set2, get prob column
################################

source_nonsmoke_Heatmap <- pathway_nonsmoke_CHL_v2 %>%
  select(source, pathway_name, prob) %>%
  group_by(source, pathway_name) %>%
  summarise(sum_total = sum(prob))

source_smoke_Heatmap <- pathway_smoke_CHL_v2 %>%
  select(source, pathway_name, prob) %>%
  group_by(source, pathway_name) %>%
  summarise(sum_total = sum(prob))

target_nonsmoke_Heatmap <- pathway_nonsmoke_CHL_v2 %>%
  select(target, pathway_name, prob) %>%
  group_by(target, pathway_name) %>%
  summarise(sum_total = sum(prob))

target_smoke_Heatmap <- pathway_smoke_CHL_v2 %>%
  select(target, pathway_name, prob) %>%
  group_by(target, pathway_name) %>%
  summarise(sum_total = sum(prob))


# transform table of source
Source_nonSK_plot <- acast(source_nonsmoke_Heatmap, pathway_name~source)
Source_SK_plot <- acast(source_smoke_Heatmap, pathway_name~source)
# transform table of target
Target_nonSK_plot <- acast(target_nonsmoke_Heatmap, pathway_name~target)
Target_SK_plot <- acast(target_smoke_Heatmap, pathway_name~target)


# now we do log2
Source_nonSK_plot_log2 <- log2(Source_nonSK_plot)
Source_SK_plot_log2 <- log2(Source_SK_plot)
Target_nonSK_plot_log2 <- log2(Target_nonSK_plot)
Target_SK_plot_log2 <- log2(Target_SK_plot)

# now get row only exist in both 
idx <- intersect(rownames(Source_nonSK_plot_log2), rownames(Source_SK_plot_log2))
idx2 <- intersect(rownames(Target_nonSK_plot_log2), rownames(Target_SK_plot_log2))

Source_nonSK_plot_log2 <- Source_nonSK_plot_log2[idx,]
Source_SK_plot_log2 <- Source_SK_plot_log2[idx,]
Target_nonSK_plot_log2 <- Target_nonSK_plot_log2[idx2,]
Target_SK_plot_log2 <- Target_SK_plot_log2[idx2,]

# log2 fold change is logp1 - logp2

Heat_plot_1 <- (Source_nonSK_plot_log2 - Source_SK_plot_log2)
Heat_plot_2 <- (Target_nonSK_plot_log2 - Target_SK_plot_log2 )

Heat_plot_1[is.na(Heat_plot_1)] <- 0
Heat_plot_2[is.na(Heat_plot_2)] <- 0

# select gene of interest
Heat_plot_1 <- Heat_plot_1[GeneInterest, , drop = FALSE]
Heat_plot_2 <- Heat_plot_2[GeneInterest, , drop = FALSE]

Heat_plot_1 <- Heat_plot_1[,cell_order]
Heat_plot_2 <- Heat_plot_2[,cell_order]

# plot heatmap

pheatmap(Heat_plot_1,cluster_rows = F,cluster_cols = F)
pheatmap(Heat_plot_2,cluster_rows = F,cluster_cols = F)

#########################################################
## for set3, first make df showing same target and source 
df_non <- pathway_nonsmoke_CHL_v2[which(pathway_nonsmoke_CHL_v2$source == pathway_nonsmoke_CHL_v2$target),]
df_non <- df_non %>%
  select(source, pathway_name, prob) %>%
  group_by(source, pathway_name) %>%
  summarise(sum_total = sum(prob))

df_non <- acast(df_non, pathway_name~source)


df_sk <- pathway_smoke_CHL_v2[which(pathway_smoke_CHL_v2$source == pathway_smoke_CHL_v2$target),]
df_sk <- df_sk %>%
  select(source, pathway_name, prob) %>%
  group_by(source, pathway_name) %>%
  summarise(sum_total = sum(prob))
df_sk <- acast(df_sk, pathway_name~source)

df_sk[is.na(df_sk)] <- 0
df_non[is.na(df_non)] <- 0

## now we me will have source mx1 and rarget mx2 before doing log2
# make NA to 0
Source_nonSK_plot[is.na(Source_nonSK_plot)] <- 0
Target_nonSK_plot[is.na(Target_nonSK_plot)] <- 0
Source_SK_plot[is.na(Source_SK_plot)] <- 0
Target_SK_plot[is.na(Target_SK_plot)] <- 0

# in Non-smoking: mx1 + mx2 - df_non, and then do log2(dataset_non)

cb_NonSk <- Source_nonSK_plot + Target_nonSK_plot 
# find interset of df_non
id_f <- intersect(rownames(cb_NonSk), rownames(df_non))
cb_NonSk <- cb_NonSk[id_f,]
df_non <- df_non[id_f,]
# no minus df_non
cb_NonSk_Final <- cb_NonSk - df_non


# in smoking: mx1 + mx2 - df_sk, and then do log2(dataset_SK)
cb_SK <- Source_SK_plot + Target_SK_plot
# find interset of df_sk
id_f2 <- intersect(rownames(cb_SK), rownames(df_sk))
cb_SK <- cb_SK[id_f2,]
df_sk <- df_sk[id_f2,]
# no minus df_non
cb_Sk_Final <- cb_SK - df_sk

# Now do lop2p
# now we select gene of interest and  do log2
# select gene of interest

cb_Sk_Final <- cb_Sk_Final[GeneInterest, , drop = FALSE]
cb_NonSk_Final <- cb_NonSk_Final[GeneInterest, , drop = FALSE]
# cell order check
cb_Sk_Final <- cb_Sk_Final[,cell_order]
cb_NonSk_Final <- cb_NonSk_Final[,cell_order]
# log2
cb_Sk_Final <- log2(cb_Sk_Final)
cb_NonSk_Final <- log2(cb_NonSk_Final)

# final,  log2(dataset_SK) -  log2(dataset_non) and then plot heatmap
set3 <-  cb_NonSk_Final - cb_Sk_Final
set3[is.na(set3)] <- NA
set3[is.infinite(as.matrix(set3))] <- NA

pheatmap(set3,cluster_rows = F,cluster_cols = F)

library(ggplot2)

set3_long <- as.data.frame(set3) %>%
  mutate(gene = rownames(set3)) %>%
  tidyr::pivot_longer(-gene, names_to = "celltype", values_to = "value")

# Convert to factors with specified order
set3_long$gene <- factor(set3_long$gene, levels = rev(GeneInterest))
set3_long$celltype <- factor(set3_long$celltype, levels = cell_order)

ggplot(set3_long, aes(x = celltype, y = gene, fill = value)) +
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", na.value = "grey80") +
  coord_fixed() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# Now do bar plot of the second figure
# Use dataset
# rankdf$relative_sum <- rankdf$contribution.scaled/max(rankdf$contribution.scaled)
# select gene of interest
# rankdf_bar <- rankdf[which(rankdf$name %in% GeneInterest),]
# rankdf_bar$name <- factor(rankdf_bar$name , levels = rev(GeneInterest))
# 
# ggplot(rankdf_bar, aes(x=relative_sum, y=name, fill=group)) +
#   geom_bar(stat='identity', position='dodge',width = 0.5, color = "black",linewidth = 0.2) +
#   scale_fill_manual(values = c("s" = "#d5d2f3", "ns" = "#ff4c4c")) + theme_classic()


# save table and plot using excel 
#write.table(rankdf_bar,'Cellchat_result/ranked_for_barplot.csv',col.names = T,row.names = F,sep = ',',quote = F)


######### Plot cellchat tutoriol plots ##########


gg1 <- netVisual_heatmap(cellchat_result_list,cluster.rows = T,cluster.cols = T)

# Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat_result_list, measure = "weight",cluster.rows = T,cluster.cols = T)
# Do heatmap based on a merged object
gg1 + gg2


# plot MHC-I and MHC-II


load('Cellchat_result/cellchat_object.list_ns_s.RData')
pathways.show <- c("MHC-I") 
par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object.list[[i]], signaling = pathways.show, color.heatmap = "Reds",title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
}

#> Do heatmap based on a single object 
#> 
#> Do heatmap based on a single object
ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))






```
