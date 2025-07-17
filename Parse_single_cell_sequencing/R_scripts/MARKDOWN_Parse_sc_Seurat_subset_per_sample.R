## Seurat 5 for parse downstream analysis in Biowulf environment ##
# date: 12052024
# Goal: making Markdown


library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(sctransform)
library(BPCells)
library(celldex)
library(SingleR)
library(biomaRt)


setwd('/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/')


######################################################################
  ############# Further analysis for diff etc... ################
######################################################################

# loading Data, for data preprocessing, find Parse_sc_Seurat_subset_per_sample.R
bladder <- readRDS('BP_FINAL_MERGEed_1M_parse_withcellanno_Intered.rds')

# check 
bladder


# plot cluster 
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
VlnPlot(object = bladder, features = 'CHRNA5',group.by = 'CELL_LAB',raster = F)
VlnPlot(bladder,features = c('CD14'),group.by = "CELL_LAB",raster = F)
VlnPlot(bladder,features = c('CD16'),group.by = "CELL_LAB",raster = F)
VlnPlot(bladder,features = c('CD68','CSF1R','FCGR3A'),group.by = "CELL_LAB",raster = F)

# Find cluster gene expression marker
# find markers for every cluster compared to all remaining cells, report only the positive ones
# check Assay
DefaultAssay(bladder)

# set idents 
Idents(bladder) <- as.factor(bladder$CELL_LAB)

gene.markers <- FindAllMarkers(bladder, only.pos = TRUE)

Top <- gene.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

VlnPlot(ref,features = "ENSG00000255260",raster = F)


#### Load TP as ref to see 

ref <- readRDS('../../../bladder_normal_from_Tabula_seruat_obj.rds')

# plot cluster 
DimPlot(ref,reduction = 'umap',label = T,raster = F,group.by = 'cell_type')
Idents(ref) <- as.factor(ref$cell_type)

# check gene marker

Ref.gene.markers <- FindAllMarkers(ref, only.pos = TRUE)

Ref.Top <- Ref.gene.markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1)

# change ENSG to geneID

# Set up biomaRt
mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl",host ="https://www.ensembl.org" )
# Get the conversion table
genes <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"),
               filters="ensembl_gene_id",
               values=Ref.Top$gene,  # Assuming Logcount has Ensembl IDs
               mart=mart)
# Create a named vector for easy conversion
gene_conv <- setNames(genes$hgnc_symbol, genes$ensembl_gene_id)
