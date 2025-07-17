
## Seurat 5 for parse downstream analysis in Biowulf environment ##
# date: 09112024
# GOAL: Deal with 1M parse single cell dataset
# sample: 75 normal bladder urinary tissue. the files are loaded from the all-sample folder
# method: using BPcell with Sketched method

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(sctransform)
library(BPCells)
library(biomaRt)
library(SingleR)


################
#Reading data###
################
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1024^8 * 1000)

###############
# Deal with Mega data from Parse
# use BPcells package to storge data on the disk, reduce pressure on memory in R
# ref: https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette.html
# https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette
###############
setwd('/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/DGE_filtered/')
# 
# # parse.data <- open_matrix_anndata_hdf5("anndata.h5ad")
# # write_matrix_dir(mat = parse.data, dir = "./parse_Mega_75_sample")
# #
# # parse.mat <- open_matrix_dir(dir = "./parse_Mega_75_sample")
# # # Read in cell meta data
# # cell_meta <- read.csv("cell_metadata.csv", row.names = 1)
# # # create obj
# # parse.object <- CreateSeuratObject(counts = parse.mat,min.features= 500, min.cells = 30, meta.data = cell_meta)
# # # SaveRDS for further use
# # saveRDS(object = parse.object,file = "./Parse_1M_75_bladder_tissue_sample.rds")
# 
# ###############
# 
# bladder <- readRDS('Parse_1M_75_bladder_tissue_sample.rds')
# # 
# # # ##########
# # # ### QC ###
# # # ##########
# # #
# bladder[["percent.mt"]] <- PercentageFeatureSet(bladder, pattern = "^MT-")
# 
# # # ##### doing QC #######
# # #
# bladder@meta.data$orig.ident <- factor("All_baldder_tissue")
# Idents(bladder) <- bladder@meta.data$orig.ident
# # #
# # # VlnPlot(bladder, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
# # # plot1 <- FeatureScatter(bladder, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # # plot2 <- FeatureScatter(bladder, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# # # plot1 + plot2
# # 
# # #subset the data
# bladder <- subset(bladder, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA < 20000 & percent.mt < 2.5)
# # 
# gc()
# # 
# bladder <- NormalizeData(bladder)
# bladder <- FindVariableFeatures(bladder, selection.method = "vst", nfeatures = 2000)


##########
##SKETCH #
##########


# bladder <- SketchData(
#   object = bladder,
#   ncells = 100000,
#   method = "LeverageScore",
#   sketched.assay = "sketch"
# )
# 
# 
# DefaultAssay(bladder) <- "sketch"
# bladder <- FindVariableFeatures(bladder, verbose = F)
# bladder <- ScaleData(bladder, verbose = F)
# bladder <- RunPCA(bladder, verbose = F)
# bladder <- FindNeighbors(bladder, dims = 1:30)
# # res change from 0.4 to 1 
# bladder <- FindClusters(bladder, resolution = 1)
# bladder <- RunUMAP(bladder, dims = 1:30, return.model = T)
# 
# print('done Sketch UMAP')

# ############################################################
# ##################### Test for SingleR  ####################
# ##################### and do Non-intergation ###############
# ############################################################
# 
# # Now read the RDS to test singleR 
# # the Parse_1M_Sketched.Rds has harmony reduction, now try to add pca reduction
# 
# # test the plot
# DimPlot(bladder, label = T, label.size = 3, reduction = "umap") + NoLegend()
# 
# # bladder[["sketch"]] <- JoinLayers(bladder[["sketch"]])
# 
# cells_in_sketch <- colnames(bladder[["sketch"]])
# 
# # bladder_sub_for_SR <- subset(bladder, cells = cells_in_sketch)
# Targett <- as.SingleCellExperiment(subset(bladder, cells = cells_in_sketch), assay = "sketch")
# 
# BRef <- readRDS('../../../bladder_normal_from_Tabula_seruat_obj.rds')
# DefaultAssay(BRef) <- "RNA"
# Logcount <- as.SingleCellExperiment(BRef, assay = "RNA")
# 
# 
# # Set up biomaRt
# mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
# 
# # Get the conversion table
# genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                filters="ensembl_gene_id", 
#                values=rownames(Logcount),  # Assuming Logcount has Ensembl IDs
#                mart=mart)
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
# 
# # Assign the new rownames
# rownames(Logcount) <- new_rownames
# # also change the GeneID from Targett, the parse gave this xxx-hg38-v130
# rownames(Targett) <- gsub("-hg38-v130","",rownames(Targett))
# 
# print('start pred')
# 
# # prediting 
# pred.grun <- SingleR(test=Targett, ref=Logcount, labels=Logcount$cell_type, de.method="wilcox")
# 
# print('done pred!')
# 
# # and put those labels back on our Seurat object and plot our on our umap.
# # after predited the cell type, can merge back to full dataset

# cell_type_predictions <- setNames(pred.grun$labels, colnames(Targett))
# # set names according to cell id
# bladder$CELL_LAB <- NA
# bladder$CELL_LAB[names(cell_type_predictions)] <- cell_type_predictions
# 
# lbls.keep <- table(pred.grun$pruned.labels)>500
# bladder$CELL_LAB <- ifelse(lbls.keep[bladder$CELL_LAB], bladder$CELL_LAB, 'Other')
# # split layer for full intergation
# # bladder[["sketch"]] <- split(bladder[["sketch"]], f = bladder$sample)
# ########################################################################################
# 
# print('start full join')
# 
# # Above code will generate xxx.full, use this as reduction input
# bladder <- ProjectData(
#   object = bladder,
#   assay = "RNA",
#   full.reduction = "pca.full",
#   sketched.assay = "sketch",
#   sketched.reduction = "pca",
#   umap.model = "umap",
#   dims = 1:30,
#   refdata = list(cellanno.full = "CELL_LAB")
# )
# 
# print('done full, start save')
# saveRDS(bladder,'Parse_1M_sketed_to_Full_cellanno.Rds')
# print('done SAVE')
# 

# plot 
# p1 <- DimPlot(bladder, reduction = "UMAP_FULL", group.by = "sample", alpha = 0.1,raster = F) + NoLegend()
# LabelClusters(plot = p1, id = "sample") 

# switch assay for other analysis using:  
#DefaultAssay(bladder) <- "RNA"
#DimPlot(bladder, label = T, label.size = 3, reduction = "UMAP_FULL", group.by = "sample", alpha = 0.1) + NoLegend()
#DimPlot(bladder, label = T, label.size = 3, reduction = "UMAP_FULL", alpha = 0.1) + NoLegend()

#############################################
# Test Pseudobulking analysis
#############################################
# Set identity classes to an existing column in meta data
# Idents(object = bladder) <- "seurat_annotations"
# table(Idents(bladder))
# bulk <- AggregateExpression(bladder, return.seurat = T, slot = "counts", assays = "RNA", group.by = c("celltype.full",
#                                                                                                      "sample", "seurat_annotations"))

# ######################################################################
#       ######### CELL ANNOTATION Using celldex ##################
# ######################################################################

# 
# #bladder[["RNA"]] <- JoinLayers(bladder[["RNA"]])
# bladder[["sketch"]] <- JoinLayers(bladder[["sketch"]])
# 
# 
# #bladder_sketch <- subset(bladder, cells = colnames(bladder[["sketch"]]))
# BD <- as.SingleCellExperiment(bladder, assay = "sketch")

#library(celldex)
#library(SingleR)
#library(dplyr)
#ref.set <- celldex::HumanPrimaryCellAtlasData()
# chcek overall
#unique(ref.set$label.main)
# check the label 
#unique(ref.set$label.fine) %>% grep('bladder',.,value = T)

# Now we have to create the SingleCellExperiment from Seurat obj
# Nomrlaized obj 
# bladder <- readRDS('DGE_filtered/saved_RDS_files/Ps_comb_lib1_nomralized.Rds')
#print('about to sce')

#sce <- as.SingleCellExperiment(bladder)

#print('done sce, next run SingleR')

# label cell 

#pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.main)

#print('done SingleR')

#pred.cnts <- SingleR::SingleR(test = BD, ref = ref.set, labels = ref.set$label.fine)
# Keep any types that have more than 10 cells to the label,
# and put those labels back on our Seurat object and plot our on our umap.
#lbls.keep <- table(pred.cnts$labels)>1000
# lbls.keep <- table(pred.cnts$labels)>10
#bladder$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')

#p1 <-DimPlot(bladder, reduction='umap', group.by='SingleR.labels',raster=T)

# ######################################################################
# ######################################################################
# Test Intergation
bladder <- readRDS('Parse_1M_sketed_to_Full_cellanno.Rds')

#p1<- DimPlot(bladder, reduction='full.umap', group.by='cellanno.full',raster=F) 

#p2 <- DimPlot(bladder, reduction='full.umap', group.by='sample',raster=F,label = T) + NoLegend() 

#p1 + p2 

#p3 <- DimPlot(bladder, 
#              reduction = "full.umap",
#              cells.highlight = WhichCells(bladder, expression = `CHRNA5` > 1),
#              cols = c('#eeeeee'),cols.highlight = c("#3f50cc"),
#              sizes.highlight = 0.2,order = T,raster = F) + NoLegend()

# original
#FeaturePlot(bladder, features = "CHRNA5-hg38-v130",order = T,raster = F)


#p4 <- FeaturePlot(bladder, features = "CHRNA5",order = T,min.cutoff =1,
#                  cells =WhichCells(bladder, expression = `CHRNA5` > 0),
#                  cols = c( "#f5cf60", "#fb3636"),raster = F,alpha = 0.8)

bladder[["sketch"]] <- split(bladder[["sketch"]], f = bladder$sample)
bladder[["RNA"]] <- split(bladder[["RNA"]], f = bladder$sample)

# check if we still use sketch 
DefaultAssay(bladder)

bladder <- FindVariableFeatures(bladder, verbose = F)
bladder <- ScaleData(bladder, verbose = F)
bladder <- RunPCA(bladder, verbose = T)
# integrate the datasets


#bladder <- IntegrateLayers(bladder, method = RPCAIntegration, orig = "pca", new.reduction = "integrated.rpca")

bladder <- IntegrateLayers(object = bladder, method = HarmonyIntegration,orig.reduction = "pca", new.reduction = "harmony",verbose = F)
bladder <- FindNeighbors(bladder, reduction = "harmony", dims = 1:30)
bladder <- FindClusters(bladder, resolution = 1, cluster.name = "harmony_clusters")
bladder <- RunUMAP(bladder, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")

print('done harmed, saving')
saveRDS(bladder,'Parse_1M_sketed_to_Full_cellanno_Harmed.Rds')
print('saved')

################################################
################################################

bladder <- readRDS('Parse_1M_sketed_to_Full_cellanno_Harmed.Rds')
bladder

# Intergation to FULL 
bladder <- ProjectIntegration(object = bladder, sketched.assay = "sketch", assay = "RNA", reduction = "harmony")
# chooe ref point 
bladder <- ProjectData(object = bladder, sketched.assay = "sketch", assay = "RNA", sketched.reduction = "harmony.full",
                      full.reduction = "harmony.full", dims = 1:30, refdata = list(celltype.full.Harmony = "cellanno.full"))


# rejoin the layrer
bladder <- RunUMAP(bladder, reduction = "harmony.full", dims = 1:30, reduction.name = "umap.Harm.full",
                  reduction.key = "UMAP_Harm_full_")
bladder <- JoinLayers(bladder)
bladder

# switch back to RNA 
DefaultAssay(bladder) <- "RNA"
p1 <- DimPlot(bladder, reduction = "umap.Harm.full", group.by = "sample", alpha = 0.1,order = T,raster = F)
p2 <- DimPlot(bladder, reduction = "umap.Harm.full", group.by = "celltype.full.Harmony", alpha = 0.1,raster = F)
p1 + p2 # + plot_layout(ncol = 1)

DimPlot(bladder, reduction = "full.umap", group.by = "sample", alpha = 0.1,order = T,raster = F)
DimPlot(bladder, reduction = "full.umap", group.by = "celltype.full.Harmony", alpha = 0.1,order = T,raster = F)

# CHRNA5 and other gene analysis
rownames(bladder) <- gsub("-hg38-v130","",rownames(bladder))

# check CHRNA5 expression per cell anno

VlnPlot(bladder,features = c('CD14','CD163','CD68','CSF1R','FCGR3A'),group.by = "celltype.full.Harmony",raster = F)


saveRDS(bladder,'Parse_1M_sketed_to_Full_cellanno_Harmed_11192024.Rds')

c10_markers <- FindMarkers(object = object, ident.1 = 10, max.cells.per.ident = 500, only.pos = TRUE)
head(c10_markers)





