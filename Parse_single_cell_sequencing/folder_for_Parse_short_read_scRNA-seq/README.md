## Making Markdown for scRNA-seq

The markdown files and scripts are currnetly located in the Biowulf data folder `parse_single_cell/R_scripts`

- parse install
- how to run
- the output
- Seurat
- cell annotation
- public database


-------

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

