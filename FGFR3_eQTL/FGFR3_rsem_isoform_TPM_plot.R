# GAOL: PLOT to see of the in house bladder tissue RNA-seq isoform TPM expression is similiar to GTEx 

library(dplyr)
library(tximport)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(xlsx)
library(edgeR)
library(ggplot2)

setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/RNA-seq/RNAseq_Normal_Urothelial_05152025/Noraml_rsem/')

edb <- EnsDb.Hsapiens.v86

tx2gene <- AnnotationDbi::select(
  edb,
  keys(edb, keytype="TXNAME"),
  columns = c("TXNAME", "GENEID"),
  keytype = "TXNAME"
)  

##################################################################################################################
### The option below is for getting isofom raw count #############################################################

names(tx2gene) <- c('transcript_id','gene_id','TXID')

lff <- list.files('.',pattern = "isoforms")
# like SetName()
names(lff) <- gsub(".isoforms.results", "", basename(lff))
# check lff to see name match with quant.sf or not

# output summary of isofrom TPM into gene TPM

txi <- tximport(
  lff,
  type    = "rsem",
  txIn = T,
  txOut = T)

# get the rawcounts
expected_counts <- txi$counts
# get the TPMs
# expected_counts <- txi$abundance

# 1. Extract the unique Ensembl gene IDs you have
gene_ids <- rownames(expected_counts)
gene_ids <- gsub('\\.[0-9]+','',gene_ids)


gene_map <- ensembldb::select(
  edb,
  keys    = gene_ids,
  keytype = "TXNAME",
  columns = c("TXNAME", "GENEID", "GENENAME")
)

# 3. Turn that into a named vector for lookup
symb <- setNames(gene_map$GENENAME, gene_map$TXNAME)

expected_counts <- as.data.frame(expected_counts)
rownames(expected_counts) <- gene_ids 

expected_counts <- expected_counts %>% mutate(GENEID=rownames(expected_counts), .before = names(expected_counts)[1] )

# Keep original Ensembl ID row names
# Add gene symbols as a separate column
expected_counts <- expected_counts %>% mutate(GENESYMBOL=ifelse(is.na(symb[gene_ids]), gene_ids, symb[gene_ids]), .before = names(expected_counts)[1] )

##################################################################################################################

# match names with sampleID

id <- readxl::read_xls('/Volumes/ifs/DCEG/Projects/DataDelivery/Prokunina/TP0325-RS7-Urothelial-Samples-RNA-seq/TP0325-RS7_QC-SUMMARY.xls')
id <- data.frame(id)
id$Vial.Label <- gsub(' ','_',id$Vial.Label)
id$Vial.Label <- gsub("_RNA$","",id$Vial.Label)

# get the name match 

names(expected_counts) <- gsub('Sample_','',names(expected_counts))

# Make a named vector for mapping
id_map <- setNames(id$Vial.Label, id$CGR.Sample.ID)

# Replace using mapping
colnames(expected_counts) <- ifelse(
  colnames(expected_counts) %in% names(id_map),
  id_map[colnames(expected_counts)],
  colnames(expected_counts)   # keep unchanged if not in mapping
)


## PLOT to see of the isoform TPM expression is similiar to GTEx 

# Now is to transfrom the data 
expected_counts_sub <- expected_counts %>% dplyr::filter(GENESYMBOL == "FGFR3")

iso_long <- expected_counts_sub %>%
  pivot_longer(
    cols = -c(GENESYMBOL, GENEID),
    names_to = "sampleid",
    values_to = "tpm"
  )

p_bladder <- iso_long %>%
  ggplot(aes(x = GENEID, y = log2(tpm+1))) +  # Changed x to isoform
  geom_boxplot(width = 0.2, linewidth = 0.2, colour = "black", 
               fill = "#A8D5E2", alpha = 0.6,outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.1, alpha = 0.3, color = "grey30") +
  stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000") +
  # Removed facet_wrap
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "none") +
  ggtitle("Bladder (n=114)")  # Add sample size manually

