### This is for Monocle3 for testing cell trad ##
### 08272024
library(monocle3)
library(dplyr)
library(Matrix)
library(Seurat)

# test one sublib
#setwd('/data/leec20/parse_single_cell/2024_output/first_run/res_SD330607_CAGATCAC-ATGTGAAG_L001/all-sample/DGE_filtered')
pdf("Monocle3_plot.pdf", width = 6, height = 4)
setwd('/data/leec20/parse_single_cell/2024_output/mega_combine_of_first_second_run/Comb_first_run/all-sample/DGE_filtered')

# load, working with large data sets

mat_path <- '.'
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if so.
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"


# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)
# check if the row are same
cell_meta <- cell_meta[colnames(mat), ]

# Create gene metadata
gene_metadata <- data.frame(
  gene_short_name = rownames(mat),
  row.names = rownames(mat)
)

# t1 <- read.csv('all_genes.csv')

# Create a Monocle3 cell_data_set object
cds <- new_cell_data_set(
  expression_data = mat,
  cell_metadata = cell_meta,
  gene_metadata = gene_metadata
)


## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 50,verbose = T)
plot_pc_variance_explained(cds)

# Reduce dimensionality and visualize the cells
cds <- reduce_dimension(cds)
# plot_cells(cds)
# check cds@colData for ploting or add column in this dataset for plotting  
# plot_cells(cds, color_cells_by="sample",label_groups_by_cluster = F)


# Group cells into cluster 
cds <- cluster_cells(cds, resolution=1e-5)
#plot_cells(cds)
plot_cells(cds, color_cells_by="partition", group_cells_by="partition")

##########################################
#Learn the trajectory graph
##########################################

#cds <- learn_graph(cds)

#plot_cells(cds,
#           color_cells_by = "partition",
#           label_groups_by_cluster=FALSE,
#           label_leaves=FALSE,
#           label_branch_points=FALSE)




##########################################
##########################################
# Find marker genes expressed by each cluster
##########################################
##########################################


marker_test_res <- top_markers(cds, group_cells_by="partition", 
                               reference_cells=1000, cores=8)

# For example, pseudo_R2 is one such measure. We can rank markers according to pseudo_R2

top_specific_markers <- marker_test_res %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(3, pseudo_R2)

top_specific_marker_ids <- unique(top_specific_markers %>% pull(gene_id))

# plot 
plot_genes_by_group(cds,
                    top_specific_marker_ids,
                    group_cells_by="partition",
                    ordering_type="cluster_row_col",
                    max.size=3)



##########################################
##########################################
# Automated annotation with Garnett
##########################################
##########################################

# First find the top markers that each annotated cell type expresses:
# assign cell type based on partitions
colData(cds)$assigned_cell_type <- as.character(partitions(cds))
# # do diff
# assigned_type_marker_test_res <- top_markers(cds,
#                                              group_cells_by="assigned_cell_type",
#                                              reference_cells=1000,
#                                              cores=8)
# 
# # Require that markers have at least JS specificty score > 0.5 and
# # be significant in the logistic test for identifying their cell type:
# # select top 5 score for per each cluster cell
# # q-val < 0.01, specific >=0.5
# 
# garnett_markers <- assigned_type_marker_test_res %>%
#   filter(marker_test_q_value < 0.01 & specificity >= 0.5) %>%
#   group_by(cell_group) %>%
#   top_n(5, marker_score)
# # Exclude genes that are good markers for more than one cell type:
# # remove duplicate
# garnett_markers <- garnett_markers %>% 
#   group_by(gene_short_name) %>%
#   filter(n() == 1)
# 
# # generate Garrnett input 
# #generate_garnett_marker_file(assigned_type_marker_test_res, file="./marker_file.txt",remove_duplicate_genes = T)
# generate_garnett_marker_file(garnett_markers, file="./marker_file.txt")
# 
library(garnett)
library(org.Hs.eg.db)

# Now train a Garnett classifier based on your marker file like this:
colData(cds)$garnett_cluster <- clusters(cds)

Bladder_classifier <- train_cell_classifier(cds = cds,
                                         marker_file = "/data/leec20/parse_single_cell/2024_output/first_run/res_SD330607_CAGATCAC-ATGTGAAG_L001/all-sample/DGE_filtered/marker_file_bladder_CellMarker2.txt",
                                         db=org.Hs.eg.db,
                                         cds_gene_id_type = "SYMBOL",
                                         num_unknown = 50,
                                         marker_file_gene_id_type = "SYMBOL",
                                         cores=8)


# Now that we've trained a classifier Bladder_classifier, we can use it to annotate the L2 cells according to type:
cds <- classify_cells(cds, Bladder_classifier,
                      db = org.Hs.eg.db,
                      cluster_extend = TRUE,
                      cds_gene_id_type = "SYMBOL")

# plot
plot_cells(cds,
           group_cells_by="partition",
           color_cells_by="cluster_ext_type")

dev.off()

print('all saveed!')
