
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
# expected_counts <- txi$counts
# get the TPMs
expected_counts <- txi$abundance

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

# Load the manifest data 
fdf <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_isoform_TMM_INT.snp.csv')
fdf <- fdf[,c(1,2)]

expected_counts_sub <- t(expected_counts_sub) %>% as.data.frame()
expected_counts_sub <- expected_counts_sub[-c(1,2),]
expected_counts_sub <- expected_counts_sub %>% mutate(Sample_Name=rownames(expected_counts_sub),.before = 1)
expected_counts_sub <- left_join(expected_counts_sub,fdf,by = 'Sample_Name')
expected_counts_sub <- expected_counts_sub[!is.na(expected_counts_sub$Predicted_Sex),]
expected_counts_sub <- expected_counts_sub[,c(11,1:10)]


# Convert from Wide to Long format
iso_long <- expected_counts_sub %>%
  pivot_longer(
    cols = starts_with("ENST"),     # Select all columns that start with ENST
    names_to = "isoform",           # Move column names to a new "isoform" column
    values_to = "tpm"            # Move the numeric values to a "counts" column
  )

iso_long$tpm <- as.numeric(iso_long$tpm)
iso_long$Predicted_Sex <- factor(iso_long$Predicted_Sex,levels = c("M","F"))

iso_long %>%
  ggplot(aes(x = isoform, y = log2(tpm + 1), fill = Predicted_Sex)) +  # 1. Added fill here
  geom_boxplot(width = 0.6, 
               linewidth = 0.5, 
               colour = "black", 
               alpha = 0.6, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) + # 2. Optional: adjust spacing
  geom_jitter(aes(group = Predicted_Sex),                  # 3. Group dots by sex
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 0.1, 
              alpha = 0.3, 
              color = "grey30") +
  stat_summary(fun = median, 
               geom = "point", 
               shape = 20, 
               size = 2, 
               color = "#FF0000",
               position = position_dodge(width = 0.8)) + # 4. Align red dots
  scale_fill_manual(values = c("M" = "#A8D5E2", "F" = "#F2BAC9")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "top")
