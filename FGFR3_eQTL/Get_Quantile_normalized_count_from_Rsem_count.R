# Date: 12292025
# Goal: Generate quantile normalized rsem expect count

library(dplyr)
library(tximport)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(xlsx)
library(edgeR)


setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/RNA-seq/RNAseq_Normal_Urothelial_05152025/Noraml_rsem/')

edb <- EnsDb.Hsapiens.v86

tx2gene <- AnnotationDbi::select(
  edb,
  keys(edb, keytype="TXNAME"),
  columns = c("TXNAME", "GENEID"),
  keytype = "TXNAME"
)  

names(tx2gene) <- c('transcript_id','gene_id','TXID')

lff <- list.files('.',pattern = ".genes.results")
# like SetName()
names(lff) <- gsub(".genes.results", "", basename(lff))
# check lff to see name match with quant.sf or not

# output summary of isofrom TPM into gene TPM

txi <- tximport(
  lff,
  type    = "rsem",
  txIn = F,
  txOut = F)

expected_counts <- txi$counts

# 1. Extract the unique Ensembl gene IDs you have
gene_ids <- rownames(expected_counts)
gene_ids <- gsub('\\.[0-9]+','',gene_ids)


gene_map <- ensembldb::select(
  edb,
  keys    = gene_ids,
  keytype = "GENEID",
  columns = c("GENEID","GENENAME")
)


# 3. Turn that into a named vector for lookup
symb <- setNames(gene_map$GENENAME, gene_map$GENEID)

expected_counts <- as.data.frame(expected_counts)
rownames(expected_counts) <- gene_ids 

expected_counts <- expected_counts %>% mutate(GENEID=rownames(expected_counts), .before = names(expected_counts)[1] )

# Keep original Ensembl ID row names
# Add gene symbols as a separate column
expected_counts <- expected_counts %>% mutate(GENESYMBOL=ifelse(is.na(symb[gene_ids]), gene_ids, symb[gene_ids]), .before = names(expected_counts)[1] )


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


# Now we save table
# write.csv(expected_counts,'../Sample_Expected_gene_counts.csv',row.names = F,quote = F)

library(edgeR)

# sample columns = everything except the two annotation columns
sample_cols <- setdiff(colnames(expected_counts), c("GENESYMBOL", "GENEID"))

# counts matrix (keep as numeric/double; do NOT force integer for RSEM expected counts)
cts <- as.matrix(expected_counts[, sample_cols])
storage.mode(cts) <- "double"
rownames(cts) <- expected_counts$GENESYMBOL

# optional filter like your old keep step
keep <- rowSums(cts) >= 10
cts <- cts[keep, , drop = FALSE]

# TMM factors
dge <- DGEList(counts = cts)
dge <- calcNormFactors(dge, method = "TMM")

# TMM-normalized log2 CPM (continuous scale)
norm_counts <- cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)

######
qnorm.DeseqNC.set <- t(apply(norm_counts, 1, rank, ties.method = "average"))
qnorm.DeseqNC.set <- qnorm(qnorm.DeseqNC.set / (ncol(qnorm.DeseqNC.set) + 1))

qnorm.DeseqNC.set <- as.data.frame(qnorm.DeseqNC.set)

# add gene annotation back (NO hard-coded 349/348)
qnorm.DeseqNC.set$GENEID <- rownames(qnorm.DeseqNC.set)
qnorm.DeseqNC.set$GENESYMBOL <- expected_counts$GENESYMBOL[match(qnorm.DeseqNC.set$GENEID,
                                                                 expected_counts$GENEID)]

# order columns: gene info first, then samples
qnorm.DeseqNC.set <- qnorm.DeseqNC.set[, c("GENESYMBOL", "GENEID", colnames(cts))]

