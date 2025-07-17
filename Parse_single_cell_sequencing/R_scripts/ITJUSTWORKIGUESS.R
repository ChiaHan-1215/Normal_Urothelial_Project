##############
# Trying to make h5ad file from parse long read data.
# Not sure this is the right way
#
##############

library(Seurat)
library(SeuratDisk)
library(SeuratObject)
library(BPCells)
library(Azimuth)
library(zellkonverter)
library(SingleCellExperiment)
library(reticulate)
library(anndata)
library(dplyr)


options(future.globals.maxSize = 1024^3 * 1000)

setwd('Long_read_pacbio_with_parse_test/ALL_32SMRT_cell_pacbio_files/aligned_BAM/NUCLEAR_PER_SAMPLE/Test_SWARM_field/WORD/')

####################

# 2. Read data
data <- open_matrix_anndata_hdf5("TESTP22.h5ad")

# write BPcell 
write_matrix_dir(
  mat = data,
  dir = "TESTP22_BP"
)

# open BPcell
mat <- open_matrix_dir(dir = "TESTP22_BP")

metadata <- LoadH5ADobs(path = "TESTP22.h5ad")
h5ad_data <- read_h5ad("TESTP22.h5ad")

# 3. Fix gene names - handle duplicates
# Convert factors to characters first!
gene_names <- as.character(h5ad_data$var$annot_gene_name)
gene_ids <- as.character(h5ad_data$var$annot_gene_id)

# Make unique names
gene_names_unique <- gene_names
is_duplicate <- duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE)
gene_names_unique[is_duplicate] <- paste0(gene_names[is_duplicate], "_", gene_ids[is_duplicate])

# 4. Create Seurat object
rownames(mat) <- gene_names_unique
seurat_obj <- CreateSeuratObject(
  counts = mat,
  meta.data = metadata,
  min.features = 200,
  min.cells = 30
)

# Verify
head(rownames(seurat_obj))
dim(seurat_obj)

#############################


sce <- readH5AD("TESTP22.h5ad")
seurat_obj <- as.Seurat(sce, counts = "X", data = "X")

# 1. Get the necessary data
counts_matrix <- GetAssayData(seurat_obj, slot = "counts")
gene_names <- seurat_obj@assays[["originalexp"]]@meta.features$annot_gene_name
metadata <- seurat_obj@meta.data  # Save the metadata

# 2. Fix gene names in counts matrix
rownames(counts_matrix) <- gene_names

# 3. Create new Seurat object with metadata
new_seurat <- CreateSeuratObject(counts = counts_matrix, 
                                 meta.data = metadata)  # Add metadata during creation

# Verify everything transferred correctly
head(new_seurat@meta.data)  # Should show dataset, sample, platform info
head(rownames(new_seurat))  # Should show gene names now
dim(new_seurat)  # Check dimensions


##############



# Read the h5ad file
data <- read_h5ad("TESTP22.h5ad")

# Get the proper gene names from var (annotations)
gene_names <- data$var$annot_gene_name

# Create count matrix with proper gene names
counts <- t(as.matrix(data$X))
rownames(counts) <- gene_names  # Set gene names before creating Seurat object

# Create Seurat object with named count matrix
data1 <- CreateSeuratObject(counts = counts, 
                            meta.data = data$obs,
                            min.features = 500, 
                            min.cells = 30)


##################

