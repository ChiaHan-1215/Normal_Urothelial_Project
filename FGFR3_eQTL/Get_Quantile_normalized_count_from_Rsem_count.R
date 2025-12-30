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
qnorm.DeseqNC.set <- qnorm.DeseqNC.set[,c(-1,-3)]

#write.table('')



### Load VcfR 
library(vcfR)
library(biomaRt)

# ==============================================================================
# 1) LOAD VCF AND EXTRACT DATA
# ==============================================================================
vcf_path <- "/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Parse_scRNA-seq/Sample_Genotyping/vcf_files/FGFR3_sub_NU_sample.vcf.gz"
vcf.file.tmp <- read.vcfR(vcf_path)

fix_df <- as.data.frame(getFIX(vcf.file.tmp))
vcf_genotypes.tmp <- extract.gt(
  vcf.file.tmp, element = "GT", return.alleles = TRUE,
  IDtoRowNames = TRUE, convertNA = TRUE
)

geno_with_pos <- cbind(
  CHR    = as.character(fix_df$CHROM),
  POS    = as.integer(fix_df$POS),
  SNP_ID = as.character(fix_df$ID),
  REF    = as.character(fix_df$REF),
  ALT    = as.character(fix_df$ALT),
  as.data.frame(vcf_genotypes.tmp, check.names = FALSE)
)

# ==============================================================================
# 2) BIOMART QUERY (hg38)
# ==============================================================================
query_chr <- gsub("chr", "", geno_with_pos$CHR, ignore.case = TRUE)
coords <- paste(query_chr, geno_with_pos$POS, geno_with_pos$POS, sep = ":")

snpMart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")

snp_results_raw <- getBM(
  attributes = c('refsnp_id', 'chrom_start', 'allele'),
  filters = 'chromosomal_region', 
  values = coords, 
  mart = snpMart
)

# Filter for exact coordinate matches
snp_results <- snp_results_raw[snp_results_raw$chrom_start %in% geno_with_pos$POS, ]

# ==============================================================================
# 3) MERGE AND STRAND DETECTION
# ==============================================================================
final_data <- merge(
  geno_with_pos, 
  snp_results, 
  by.x = "POS", 
  by.y = "chrom_start",
  all.x = TRUE
)

# Function to get the complement (A <-> T, C <-> G)
complement_seq <- function(sequence) {
  # Handles single bases or strings like "A/C"
  chartr("ATCG", "TAGC", sequence)
}

# Identify Strand Status
final_data$Strand_Status <- sapply(1:nrow(final_data), function(i) {
  vcf_ref <- final_data$REF[i]
  ens_alleles <- final_data$allele[i]
  if (is.na(ens_alleles)) return("Unknown")
  
  if (grepl(vcf_ref, ens_alleles)) {
    return("Forward")
  } else if (grepl(complement_seq(vcf_ref), ens_alleles)) {
    return("Reverse/Flipped")
  } else {
    return("Mismatch")
  }
})

# ==============================================================================
# 4) AUTOMATIC GENOTYPE UNF-FLIPPING (The Fix)
# ==============================================================================
# This section iterates through all sample columns and complements the alleles
# ONLY if the Strand_Status is "Reverse/Flipped".

sample_names <- colnames(vcf_genotypes.tmp)

for (sample in sample_names) {
  # Apply complement only where flipped
  # sample <- sample_names[2]
  final_data[[sample]] <- ifelse(
    final_data$Strand_Status == "Reverse/Flipped",
    complement_seq(final_data[[sample]]), 
    final_data[[sample]]
  )
}

# ==============================================================================
# 5) FINAL CLEANUP AND EXPORT
# ==============================================================================
# Update IDs to official rsIDs
final_data$Clean_SNP_ID <- ifelse(
  !is.na(final_data$refsnp_id), final_data$refsnp_id, final_data$SNP_ID
)

# Rename Ensembl column
colnames(final_data)[colnames(final_data) == "allele"] <- "Ensembl_Alleles_Plus_Strand"

# Organize columns: Metadata first, then Samples
info_cols <- c("CHR", "POS","SNP_ID","Clean_SNP_ID", "REF", "ALT", 
               "Ensembl_Alleles_Plus_Strand", "Strand_Status")
final_data <- final_data[, c(info_cols, sample_names)]

# Keep the final data that both REF and ALT are not NA
final_data <- final_data %>%
  dplyr::filter(!is.na(REF) & !is.na(ALT))


# modified the sample name 


colnames(qnorm.DeseqNC.set)
colnames(final_data)

# load sample manifest 
gt_sample_man <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Parse_scRNA-seq/Sample_Genotyping/SR0325-017_1_QC_Report_48202503941.xlsx')

gt_sample_man$Sample_Name <- gsub(" ","_",gt_sample_man$Sample_Name)
gt_sample_man$Sample_Name <- gsub("_DNA","",gt_sample_man$Sample_Name)

# cahnge name 

guide <- setNames(gt_sample_man$Sample_Name,gt_sample_man$Sample_ID)

old_names <- colnames(final_data)[9:ncol(final_data)]

# Replace only where a match exists in guide
new_names <- ifelse(
  old_names %in% names(guide),
  guide[old_names],
  old_names
)

colnames(final_data)[9:ncol(final_data)] <- new_names


# Now is to transfrom the data 

# Extract FGFR3 

FGFR3_qn <- qnorm.DeseqNC.set %>% dplyr::filter(GENEID == "FGFR3")
FGFR3_qn <- t(FGFR3_qn) %>% as.data.frame()
FGFR3_qn$Sample_Name <- rownames(FGFR3_qn)
FGFR3_qn <- FGFR3_qn[-1,]
FGFR3_qn <- FGFR3_qn[,c(2,1)]

# Test it use rsid only 
final_data_sub <- final_data[grep("rs",final_data$SNP_ID),]
final_data_sub <- final_data_sub[,c(4,9:128)]
final_data_sub <- t(final_data_sub) %>% as.data.frame()
names(final_data_sub) <- final_data_sub[1,]
final_data_sub <- final_data_sub[-1,]

# remove duplicateion 
final_data_sub <- final_data_sub[,!duplicated(names(final_data_sub))]

final_data_sub <- final_data_sub %>% mutate(Sample_Name=rownames(final_data_sub),.before = 1)

FGFR3_Merged_data <- left_join(FGFR3_qn,final_data_sub,by="Sample_Name")

# Load manifest file of sample age sex etc 
gt_sample_man_sub <- gt_sample_man %>% dplyr::select(Sample_Name,Predicted_Sex,Ancestry)

final_FGFR3 <- left_join(gt_sample_man_sub,FGFR3_Merged_data,by = 'Sample_Name')

# Done save 
# write.table(final_FGFR3,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Parse_scRNA-seq/Sample_Genotyping/FGFR3_qn.snp.Test.csv',row.names = F,quote = F,sep = ',')


  
