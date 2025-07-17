## Seurat 5 for parse downstream analysis in Biowulf environment ##
# date: 11192024
# Goal: To load each sample indivalually, and do QC and cell annotation per sample and then merge it

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(sctransform)
library(BPCells)
library(celldex)
library(SingleR)
library(biomaRt)


################
#Reading data###
################
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1024^3 * 1000)

# setwd('/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/')
# 
# ################################# load ref for annotation ################################################
# BRef <- readRDS('../bladder_normal_from_Tabula_seruat_obj.rds')
# DefaultAssay(BRef) <- "RNA"
# Logcount <- as.SingleCellExperiment(BRef, assay = "RNA")
# 
# # Set up biomaRt
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host ="https://www.ensembl.org" )
# # Get the conversion table
# genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                filters="ensembl_gene_id", 
#                values=rownames(Logcount),  # Assuming Logcount has Ensembl IDs
#                mart=mart)
# 
# # Create a named vector for easy conversion
# gene_conv <- setNames(genes$hgnc_symbol, genes$ensembl_gene_id)
# # Convert the rownames of Logcount, keeping original names if no match found
# new_rownames <- sapply(rownames(Logcount), function(x) {
#   if (x %in% names(gene_conv) && gene_conv[x] != "") {
#     return(gene_conv[x])
#   } else {
#     return(x)  # Keep the original name if no match or empty conversion
#   }
# })
# 
# # Ensure uniqueness of new rownames
# new_rownames <- make.unique(new_rownames)
# # Assign the new rownames
# rownames(Logcount) <- new_rownames

# ################################################################################################
# 
# setwd('/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/')
# 
# sample_dirs <- list.dirs(path = ".", full.names = TRUE, recursive = FALSE)
# sample_dirs <- grep("sample_[3-8][0-9]$", sample_dirs, value = TRUE)
# sample_dirs <- sample_dirs[-c(1,2)]
# 
# seurat_objects <- list()
# 
# for (i in sample_dirs[31:47]) {
#   # Get sample name from directory name
#   # i <- sample_dirs[1]
#   sample_name <- basename(i)
#   # Set paths for matrix and metadata
#   mat_path <- file.path(i, "DGE_filtered")
#   # Read ParseBio matrix
#   mat <- ReadParseBio(mat_path)
#   # Check and fix empty gene names
#   empty_genes <- rownames(mat) == ""
#   rownames(mat)[empty_genes] <- paste0("unknown_", seq_len(sum(empty_genes)))
#   # Read cell metadata
#   cell_meta <- read.csv(file.path(mat_path, "cell_metadata.csv"), row.names = 1)
#   # Create Seurat object
#   seurat_obj <- CreateSeuratObject(mat, min.features = 200, min.cells = 3, meta.data = cell_meta)
#   # Store Seurat object in list
#   seurat_objects[[sample_name]] <- seurat_obj
#   print('done one!!!')
# }
# 
# 
# # Merge all Seurat objects
# bladder <- merge(seurat_objects[[1]], y = seurat_objects[-1] , project = "combined_bladder")
# # join layer if merge
# bladder[["RNA"]] <- JoinLayers(bladder[["RNA"]])
# 
# 
# # ##########
# # ### QC ###
# # ##########
# #
# bladder[["percent.mt"]] <- PercentageFeatureSet(bladder, pattern = "^MT-")
# 
# # ##### doing QC #######
# #
# bladder@meta.data$orig.ident <- factor("All_baldder_tissue")
# Idents(bladder) <- bladder@meta.data$orig.ident
# #
# VlnPlot(bladder, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# 
# plot1 <- FeatureScatter(bladder, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(bladder, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# 
# # # # # subset the data 
# 
# bladder <- subset(bladder, subset = nFeature_RNA > 200 & nFeature_RNA < 12500 & nCount_RNA < 12500 & percent.mt < 35)
# gc()
# 
# bladder <- NormalizeData(bladder)
# bladder <- FindVariableFeatures(bladder, selection.method = "vst", nfeatures = 2000)
# all.genes <- rownames(bladder)
# bladder <- ScaleData(bladder, features = all.genes)
# bladder <- RunPCA(bladder, features = VariableFeatures(object = bladder))
# 
# ElbowPlot(bladder)
# 
# bladder <- FindNeighbors(bladder, dims = 1:30)
# bladder <- FindClusters(bladder,resolution = 0.5)
# bladder <- RunUMAP(bladder, dims = 1:30)
# 
# 
# # ######################################################################
# # ################################### CELL ANNOTATION ##################
# # ######################################################################

# # set martix
# Targett <- as.SingleCellExperiment(bladder, assay = "RNA")
# 
# ####### prediting
# pred.grun <- SingleR(test=Targett, ref=Logcount, labels=Logcount$cell_type, de.method="wilcox")
# 
# # and put those labels back on our Seurat object and plot our on our umap.
# 
# lbls.keep <- table(pred.grun$pruned.labels)>10
# 
# bladder$CELL_LAB <- ifelse(lbls.keep[pred.grun$labels], pred.grun$labels, 'Other')
# 
# # DimPlot(bladder,group.by = 'CELL_LAB',reduction = 'umap',raster = F,label = T) 
# print('done anno, prepate save')
# # save it as h5 
# #TarS <- bladder$sample %>% unique()
# #TarS
# 
# saveRDS(bladder,"/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/Merge_67_end.rds")
# 
# print('done save')

######################################################################
############# 10222024 Mege huge rds file ############################
######################################################################

# Set working directory and get file list
setwd('/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/')
# sample_ff <- list.files(full.names = TRUE, pattern = '.rds')
# 
# # Initialize with first object
# print(paste("Starting merge of", length(sample_ff), "files"))
# 
# bladder <- readRDS(sample_ff[1])
# 

# Remove unused field for reduce memory issue

# bladder[["RNA"]]$data <- NULL
# bladder[["RNA"]]$scale.data <- NULL
# bladder[["RNA"]]@meta.data <- data.frame(row.names = rownames(bladder[["RNA"]]))
# bladder[["pca"]] <- NULL
# bladder[["umap"]] <- NULL
# 
# gc()
# 
# # Merge remaining files
# for(i in 2:length(sample_ff)) {
#   print(paste("Merging file", i, "of", length(sample_ff)))
#   current_obj <- readRDS(sample_ff[i])
#   current_obj[["RNA"]]$data <- NULL
#   current_obj[["RNA"]]$scale.data <- NULL
#   current_obj[["RNA"]]@meta.data <- data.frame(row.names = rownames(current_obj[["RNA"]]))
#   current_obj[["pca"]] <- NULL
#   current_obj[["umap"]] <- NULL
#   
#   bladder <- merge(bladder,y = current_obj)
#   rm(current_obj); gc()
# }
# 
# bladder<-JoinLayers(bladder)
# 

# write up the BPcell storage 

# write_matrix_dir(mat = bladder[["RNA"]]$counts, dir = '/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/BPcell_storage_for_all/')
# counts.mat <- open_matrix_dir(dir = '/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/BPcell_storage_for_all/')
# bladder[["RNA"]]$counts <- counts.mat
# 

######################################################################
############# Do with or without Integration analysis ################
######################################################################

bladder <- readRDS("BP_FINAL_MERGEed_1M_parse_withcellanno.rds")

bladder[["RNA"]] <- split(bladder[["RNA"]], f = bladder$sample)

# remove unused column
columns_to_remove <- c("RNA_snn_res.0.1")
bladder@meta.data <- bladder@meta.data[, !colnames(bladder@meta.data) %in% columns_to_remove]
# 
bladder <- NormalizeData(bladder)
bladder <- FindVariableFeatures(bladder)
bladder <- ScaleData(bladder)
bladder <- RunPCA(bladder,verbose = F)

# 
## Without interg
bladder <- FindNeighbors(bladder, dims = 1:30, reduction = "pca")
# # res can set 0.1-2, depend on how you want seurat cluster number
bladder <- FindClusters(bladder, resolution = 1, cluster.name = "unintegrated_clusters")
#bladder <- FindClusters(bladder, resolution = 0.5, cluster.name = "unintegrated_clusters")
bladder <- RunUMAP(bladder, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# DimPlot(bladder, reduction = "umap.unintegrated")
# VlnPlot(bladder,features = 'CHRNA5',group.by = 'CELL_LAB')

## With interg
bladder <- IntegrateLayers(
  object = bladder, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE)

# # re-join after interger
bladder <- JoinLayers(bladder)

bladder <- FindNeighbors(bladder, reduction = "harmony", dims = 1:30)
bladder <- FindClusters(bladder, resolution = 1, cluster.name = "harmony_clusters")
bladder <- RunUMAP(bladder, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

#saveRDS(bladder, "BP_FINAL_MERGEed_1M_parse_withcellanno_Intered.rds")

######################################################################
  ############# Further analysis for diff etc... ################
######################################################################

bladder <- readRDS('BP_FINAL_MERGEed_1M_parse_withcellanno_Intered.rds')
bladder




DimPlot(bladder,reduction = 'umap.harmony',label = T,group.by = "sample",raster = F)
DimPlot(bladder,reduction = 'umap.harmony',label = T,group.by = "CELL_LAB",raster = F,order = T,alpha = 0.7) + NoLegend()

DimPlot(bladder,reduction = 'umap.unintegrated',label = T,group.by = "sample",raster = F) + NoLegend()
DimPlot(bladder,reduction = 'umap.unintegrated',label = T,group.by = "CELL_LAB",raster = F,order = T,alpha = 0.7)



FeaturePlot(bladder, features = "CHRNA5",reduction = "umap.harmony",order = T,min.cutoff =1,
                  cells =WhichCells(bladder, expression = `CHRNA5` > 0),
                  cols = c( "#f5cf60", "#fb3636"),raster = F,alpha = 0.8)


DimPlot(bladder, 
              reduction = "umap.harmony",
              cells.highlight = WhichCells(bladder, expression = `CHRNA5` > 1),
              cols = c('#eeeeee'),cols.highlight = c("#3f50cc"),
              sizes.highlight = 0.2,order = T,raster = F) + NoLegend()

# testing gene marker
###########################################################################
############          SEX detemine #######
###########################################################################

### Subset sample so that it's earsier to plot

# Get all unique sample names
all_samples <- unique(bladder$sample)
# Find those that start with "P_"
p_samples <- grep("^P", all_samples, value = TRUE)
Other_samples <- grep("sample_", all_samples, value = TRUE)
# subset 15

sex_fea <- c('RPS4Y1','EIF1AY','UTY','XIST','KDM6A')



for ( i in seq(1, 75, by=15)){
  bladder_subset <- subset(bladder, subset = sample %in% all_samples[i:min(i+14, 75)])
  DotPlot(bladder_subset,features = c('UTY','XIST'),group.by = 'sample')
  rm(bladder_subset)
}

# 
# bladder_subset <- subset(bladder, subset = sample %in% Other_samples[1:15])
# ## UTY and XIST for M and F is good marker ##
# DotPlot(bladder_subset,features = c('UTY','XIST'),group.by = 'sample')
# VlnPlot(bladder_subset,features =c('UTY','XIST'),group.by = 'sample',pt.size = 0.1,ncol = 2)



###########################################################################
###########################################################################


VlnPlot(bladder,features = c('CD14'),group.by = "CELL_LAB",raster = F)
VlnPlot(bladder,features = c('CD16'),group.by = "CELL_LAB",raster = F)
VlnPlot(bladder,features = c('CD68','CSF1R','FCGR3A'),group.by = "CELL_LAB",raster = F)

# Find cluster gene expression marker
# find markers for every cluster compared to all remaining cells, report only the positive ones


pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE)
pbmc.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)







# Finding which cell high expressed CHRAN5 vs low express, and finding coexpression gene #

#chrna5_expr <- FetchData(bladder, vars = "CHRNA5")
#bladder$CHRNA5_groups <- ifelse(chrna5_expr$CHRNA5 >= 1.5, "CHRNA5_high", "CHRNA5_low")
#Idents(bladder) <- "CHRNA5_groups"
#VlnPlot(bladder, features = c("SCN11A"), group.by = "CHRNA5_groups")

# de_results <- FindMarkers(bladder,
#                           ident.1 = "CHRNA5_high",
#                           ident.2 = "CHRNA5_low",
#                           min.pct = 0.25,
#                           test.use = "wilcox")


