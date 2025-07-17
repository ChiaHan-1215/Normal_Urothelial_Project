## Seurat 5 for parse downstream analysis in Biowulf environment ##
# date: 08152024

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(sctransform)
library(BPCells)
################
#Reading data###
################
# needs to be set for large dataset analysis
options(future.globals.maxSize = 1024^3 * 1000)

###############
# Deal with Mega data from Parse
# ref: https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette.html
# https://satijalab.org/seurat/articles/seurat5_bpcells_interaction_vignette
###############
setwd('/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/DGE_filtered/')
#setwd('/data/leec20/parse_single_cell/2024_output/mega_combine_of_first_second_run/Comb_first_run/all-sample/')
# 
# parse.data <- open_matrix_anndata_hdf5("anndata.h5ad")
# write_matrix_dir(mat = parse.data, dir = "./parse_Mega_75_sample")
# 
# parse.mat <- open_matrix_dir(dir = "./parse_Mega_75_sample")
# # Read in cell meta data
# cell_meta <- read.csv("cell_metadata.csv", row.names = 1)
# # create obj
# parse.object <- CreateSeuratObject(counts = parse.mat,min.features= 500, min.cells = 30, meta.data = cell_meta)
# # SaveRDS for further use 
# saveRDS(object = parse.object,file = "./Parse_1M_75_bladder_tissue_sample.rds")

bladder <- readRDS('Parse_1M_75_bladder_tissue_sample.rds')

# # # 
# mat_path <- "DGE_filtered/"
# mat <- ReadParseBio(mat_path)

# # Check to see if empty gene names are present, add name if so.
# table(rownames(mat) == "")
# rownames(mat)[rownames(mat) == ""] <- "unknown"
# 
# # Read in cell meta data
# cell_meta <- read.csv(paste0(mat_path, "cell_metadata.csv"), row.names = 1)
# 
# # create Sureat object
# bladder <- CreateSeuratObject(mat, min.features= 200, min.cells = 3, meta.data = cell_meta)
# Setting our initial cell class to a single type, this will changer after clustering.
#

# ##########
# ### QC ###
# ##########
#
bladder[["percent.mt"]] <- PercentageFeatureSet(bladder, pattern = "^MT-")
#
# ###### OPTIONAL: Visualize QC metrics as a violin plot ######
#
# ## Visual how many cell counts per group ######
# bladder@meta.data$orig.ident <- factor(bladder@meta.data$sample)
# Idents(bladder) <- bladder@meta.data$orig.ident
#
# # bladder@meta.data %>% ggplot(aes(x=orig.ident, fill=orig.ident)) +
# #   geom_bar(color="black") + stat_count(geom = "text", colour = "black", size = 3.5, aes(label = ..count..),position=position_stack(vjust=0.5))+
# #   theme_classic() + theme(plot.title = element_text(hjust=0.5, face="bold")) +
# #   ggtitle("Number of Cells per Sample")
#
# ##### doing QC #######
#
bladder@meta.data$orig.ident <- factor("All_baldder_tissue")
Idents(bladder) <- bladder@meta.data$orig.ident
#
# VlnPlot(bladder, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# plot1 <- FeatureScatter(bladder, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(bladder, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1 + plot2
# # 
# # # # subset the data 
bladder <- subset(bladder, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & nCount_RNA < 20000 & percent.mt < 3)
# #
gc()

#############################################
# Normalizing the data using SCtransform V2 # 
#############################################

# Around 750G is save after combine all data for traditional normalization, SCT is too heavey memory which use over 800G...
# read Rds
# bladder <- readRDS('DGE_filtered/saved_RDS_files/Ps_comb_lib1_nomralized.Rds')
# bladder <- SCTransform(bladder, vars.to.regress = "percent.mt", verbose = FALSE,ncells = 1000)
# bladder <- SCTransform(bladder, vars.to.regress = "percent.mt", verbose = FALSE)

# print('SCTed')

# # run PCA and other stuff
# bladder <- RunPCA(bladder, verbose = FALSE,assay = "SCT")
# quick check how many pcs we need
# ElbowPlot(bladder, ndims = 40)
# use Seruat recommand 30PCs
# bladder <- FindNeighbors(bladder, dims = 1:30)
# bladder <- FindClusters(bladder,resolution = 0.1)
# bladder <- RunUMAP(bladder, dims = 1:30)


# SAVE NoW 
# saveRDS(bladder,file = 'DGE_filtered/saved_RDS_files/Ps_comb_lib1_SCT_done.Rds')


########################################################################################################################################################################
########################################################################################################################################################################
# traditional way 

bladder <- NormalizeData(bladder)
bladder <- FindVariableFeatures(bladder, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(bladder)
bladder <- ScaleData(bladder, features = all.genes)

print('done Scaled')

bladder <- RunPCA(bladder, features = VariableFeatures(object = bladder))
# ElbowPlot(bladder)
bladder <- FindNeighbors(bladder, dims = 1:30)

print('Done PCA')

bladder <- FindClusters(bladder,resolution = 0.1)
bladder <- RunUMAP(bladder, dims = 1:30)

print('about to save')

saveRDS(bladder,file = 'Parse_Non_intergated_1M_75_bladder_tissue_sample.rds')
print('done saved')


#################################################################################
# Selecting cells for plotting or diff analysis  ##########################
#################################################################################
# Once we saved, now load the file, it's ~150G
#bladder <- readRDS('Parse_Non_intergated_1M_75_bladder_tissue_sample.rds')

# 
# # 
# # check assay
# #Check default assay
# DefaultAssay(object = bladder)
# # 
# # # plot
DimPlot(bladder, label = TRUE,reduction = 'umap',group.by = 'sample', alpha = 0.1,raster=FALSE) + NoLegend()
DimPlot(bladder, label = TRUE,reduction = 'umap', alpha = 0.1,raster=F) + NoLegend()
## Plot CHRNA5 expression 
rownames(bladder) <- gsub("-hg38-v130",'',rownames(bladder))

VlnPlot(bladder, features = "CHRNA5",raster = F) + NoLegend()

# Create a UMAP plot CHRNA5-expressing cells, not highlight

p1 <- DimPlot(bladder, 
        reduction = "umap",
        cells.highlight = WhichCells(bladder, expression = `CHRNA5` > 1),
        cols = c('#eeeeee'),cols.highlight = c("#3f50cc"),
        sizes.highlight = 0.2,order = T,raster = F) + NoLegend()


# original
#FeaturePlot(bladder, features = "CHRNA5-hg38-v130",order = T,raster = F)


p2 <- FeaturePlot(bladder, features = "CHRNA5",order = T,min.cutoff =1,
            cells =WhichCells(bladder, expression = `CHRNA5` > 0),
            cols = c( "#f5cf60", "#fb3636"),raster = F,alpha = 0.8)




# # or by functuion WhichCells(), by idents 
# cells_to_plot <- WhichCells(bladder, idents = c("12", "14"))
# # by sample_name etc..
# cells_P7<- WhichCells(bladder, expression = sample == "P7")

# ### Additionally, only plot the sample of cell 
# DimPlot(bladder, 
#         reduction = "umap",
#         group.by = "seurat_clusters",
#         cells = cells_P11) +
#   labs(title = "UMAP plot of Sample P11 colored by clusters")

############ For mamually select cell and do pseudo bulk analysis

# # select hightlight cells of interest by window
# pp <- DimPlot(bladder, reduction = "umap",alpha = 0.1,raster=FALSE)
# select.cells <- CellSelector(plot = pp)
# 
# # now create new group 
# bladder$selected_group <- "Other"
# bladder$selected_group[select.cells] <- "Selected"
# 
# # now test the pseudobulk
# 
# pseudobulk <- AggregateExpression(
#   bladder,
#   group.by = c("selected_group", "sample"),
#   assays = "RNA",
#   return.seurat = T
# )
# 
# # check each 'cell' is a donor-condition-celltype pseudobulk profile
# tail(Cells(pseudobulk))
# 
# Idents(pseudobulk) <- "selected_group"
# 
# bulk.mono.de <- FindMarkers(object = pseudobulk, 
#                             ident.1 = "Selected", 
#                             ident.2 = "Other",
#                             test.use = "DESeq2")

######## ######## ######## ######## ########
####### Idents ########

# set Idents
#Idents(bladder) <- bladder@meta.data$sample

# rename Idents
#bladder <- RenameIdents(bladder, 'P11' = 'SAMPLE_P11', 'P14_cult' = 'P14_culture')

# plot gene expression for testing
# Visualize canonical marker genes as violin plots.
# VlnPlot(bladder, features = c("SNORD3B-1"),idents = c("P12","P15"),pt.size = 0.2)

######## ######## ######## ######## ########
######## ######## ######## ######## ########

# ######################################################################
# ###########                         Expression Diff ##################
# ######################################################################
# 
# # find all markers of cluster 2
# 
# library(presto)
# 
# # # find all marker 
# bladder.markers <- FindAllMarkers(bladder, only.pos = TRUE,logfc.threshold = 0.5)
# 
# bladder.markers %>%
#   group_by(cluster) %>%
#   dplyr::filter(avg_log2FC > 1)
# # 
# # #ordering the results 
# bladder.markers <- bladder.markers %>% arrange(cluster,desc(avg_log2FC), desc(p_val_adj))
# # #examine a small subset 
# bladder.markers %>% group_by(cluster) %>% slice_max(n = 5, order_by = avg_log2FC)
# # 
# #plot(density(sample(JoinLayers(bladder@assays$RNA)$count["ENSG00000290457",],2500)))
# # 
# VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

### Seting Ident to compare different group
# ifnb$celltype.stim <- paste(ifnb$seurat_annotations, ifnb$stim, sep = "_")
# Idents(ifnb) <- "celltype.stim"
# mono.de <- FindMarkers(ifnb, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", verbose = FALSE)
# head(mono.de, n = 10)

# ######################################################################
# ################################### CELL ANNOTATION ##################
# ######################################################################
# using Azimuth? 
# Do cell intergation
# Try other stuff 
library(celldex)
library(SingleR)
library(dplyr)
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

#pred.cnts <- SingleR::SingleR(test = sce, ref = ref.set, labels = ref.set$label.fine)
# Keep any types that have more than 10 cells to the label,
# and put those labels back on our Seurat object and plot our on our umap.
#lbls.keep <- table(pred.cnts$labels)>1000
# lbls.keep <- table(pred.cnts$labels)>10
#bladder$SingleR.labels <- ifelse(lbls.keep[pred.cnts$labels], pred.cnts$labels, 'Other')

#p1 <-DimPlot(bladder, reduction='umap', group.by='SingleR.labels',raster=T)


# save 
#pdf("SingleR_dimplot.pdf", width = 12, height =10)
#print(p1)


#dev.off()
#print('done Save plot')
# 
# #############################################
# #############################################
# # prepare count slots for seurat obj
# # If you need to create a counts slot for compatibility with other functions
# bladder[["RNA3"]] <- as(object = bladder[["RNA"]], Class = "Assay")
# 
# # CellAnn Package
# 
# ### download scripts ###
# 
# devtools::source_url("https://raw.githubusercontent.com/Pinlyu3/CellAnn/main/prepare_CellAnn.R")
# 
# ### prepare_CellAnn function:
# ### parameter: seurat_obj: your seurat obj
# ### parameter: folder(character): CellAnn ouput files will be output to this path, default is your current folder
# ### parameter: sample_name(character): names of this sample, such as "Liver_1", "eye2_2" etc.
# ### parameter: matrix_name(character): default will use the "RNA" counts matrix in your seurat object
# ### parameter: dims(character): "umap" or "tsne". default is 'umap'
# ### parameter: cluster(character): the name of the column which stored cluster information in the metadata of your Seurat object. default is 'seurat_clusters'
# 
# prepare_CellAnn(bladder,folder='/data/leec20/parse_single_cell/2024_output/first_run/Comb_first_run/all-sample/CellAnn/',sample_name='CellAnn_cblib1',matrix_name='RNA3',dims='umap',cluster='seurat_clusters')
# 
# ### After run prepare_CellAnn function, you will find 2 prepared files under your folder 
# # need to remove ENSG stuff in our data 
# df <- read.table('CellAnn_cblib1_CellAnn_Step1_input.txt',header = T,check.names = FALSE)
# # remvoe ENSG Gene
# df <- df %>% filter(!grepl("^ENSG", GENE))
# write.table(df,'CellAnn_cblib1_CellAnn_Step1_input_md.txt',col.names = T,row.names = F,sep = '\t',quote = F)
# 
# 
# #############################################
# #############################################

BRef <- readRDS('../../../bladder_normal_from_Tabula_seruat_obj.rds')
DefaultAssay(BRef) <- "RNA"

# set martix
Targett <- as.SingleCellExperiment(bladder, assay = "RNA")
Logcount <- as.SingleCellExperiment(BRef, assay = "RNA")

######
library(biomaRt)
# Set up biomaRt
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# Get the conversion table
genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
               filters="ensembl_gene_id", 
               values=rownames(Logcount),  # Assuming Logcount has Ensembl IDs
               mart=mart)


# Create a named vector for easy conversion
gene_conv <- setNames(genes$hgnc_symbol, genes$ensembl_gene_id)

# Convert the rownames of Logcount, keeping original names if no match found
new_rownames <- sapply(rownames(Logcount), function(x) {
  if (x %in% names(gene_conv) && gene_conv[x] != "") {
    return(gene_conv[x])
  } else {
    return(x)  # Keep the original name if no match or empty conversion
  }
})


# Ensure uniqueness of new rownames
new_rownames <- make.unique(new_rownames)

# Assign the new rownames
rownames(Logcount) <- new_rownames

####### prediting

pred.grun <- SingleR(test=Targett, ref=Logcount, labels=Logcount$cell_type, de.method="wilcox")

# and put those labels back on our Seurat object and plot our on our umap.

lbls.keep <- table(pred.grun$pruned.labels)>500

bladder$CELL_LAB <- ifelse(lbls.keep[pred.grun$labels], pred.grun$labels, 'Other')


DimPlot(BRef,group.by = 'cell_type',reduction = 'umap')
DimPlot(bladder,group.by = 'CELL_LAB',reduction = 'umap',raster = F)

