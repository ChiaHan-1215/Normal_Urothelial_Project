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

# lff <- list.files('.',pattern = ".genes.results")
# # like SetName()
# names(lff) <- gsub(".genes.results", "", basename(lff))
# # check lff to see name match with quant.sf or not
# 
# # output summary of isofrom TPM into gene TPM
# 
# txi <- tximport(
#   lff,
#   type    = "rsem",
#   txIn = F,
#   txOut = F)
# 
# expected_counts <- txi$counts
# 
# # 1. Extract the unique Ensembl gene IDs you have
# gene_ids <- rownames(expected_counts)
# gene_ids <- gsub('\\.[0-9]+','',gene_ids)
# 
# 
# gene_map <- ensembldb::select(
#   edb,
#   keys    = gene_ids,
#   keytype = "GENEID",
#   columns = c("GENEID","GENENAME")
# )
# 
# 
# # 3. Turn that into a named vector for lookup
# symb <- setNames(gene_map$GENENAME, gene_map$GENEID)
# 
# expected_counts <- as.data.frame(expected_counts)
# rownames(expected_counts) <- gene_ids 
# 
# expected_counts <- expected_counts %>% mutate(GENEID=rownames(expected_counts), .before = names(expected_counts)[1] )
# 
# # Keep original Ensembl ID row names
# # Add gene symbols as a separate column
# expected_counts <- expected_counts %>% mutate(GENESYMBOL=ifelse(is.na(symb[gene_ids]), gene_ids, symb[gene_ids]), .before = names(expected_counts)[1] )
# 

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


# Now we save table
# write.csv(expected_counts,'../Sample_Expected_gene_counts.csv',row.names = F,quote = F)

library(edgeR)

# sample columns = everything except the two annotation columns
sample_cols <- setdiff(colnames(expected_counts), c("GENESYMBOL", "GENEID"))

# counts matrix (keep as numeric/double; do NOT force integer for RSEM expected counts)
cts <- as.matrix(expected_counts[, sample_cols])
storage.mode(cts) <- "double"
rownames(cts) <- expected_counts$GENEID

# optional filter like your old keep step
keep <- rowSums(cts) >= 10
cts <- cts[keep, , drop = FALSE]

###############################################################################################################
# Option1: TMM then do quantile normalized
###############################################################################################################

# # TMM factors
# dge <- DGEList(counts = cts)
# dge <- calcNormFactors(dge, method = "TMM")
# # TMM-normalized log2 CPM (continuous scale)
# norm_counts <- cpm(dge, log = TRUE, prior.count = 1, normalized.lib.sizes = TRUE)
# 
# ######
# qnorm.DeseqNC.set <- t(apply(norm_counts, 1, rank, ties.method = "average"))
# qnorm.DeseqNC.set <- qnorm(qnorm.DeseqNC.set / (ncol(qnorm.DeseqNC.set) + 1))
# 
# qnorm.DeseqNC.set <- as.data.frame(qnorm.DeseqNC.set)
# 
# # add gene annotation back (NO hard-coded 349/348)
# qnorm.DeseqNC.set$GENEID <- rownames(qnorm.DeseqNC.set)
# qnorm.DeseqNC.set$GENESYMBOL <- expected_counts$GENESYMBOL[match(qnorm.DeseqNC.set$GENEID,
#                                                                  expected_counts$GENEID)]
# 
# # order columns: gene info first, then samples
# qnorm.DeseqNC.set <- qnorm.DeseqNC.set[, c("GENESYMBOL", "GENEID", colnames(cts))]
# qnorm.DeseqNC.set <- qnorm.DeseqNC.set[,c(-1,-3)]

###############################################################################################################
# Option2: To follow GTEx v8-v10 version, they use TMM then do rank-based Inverse Normal Transformation (INT)
# Extract from Oscar's script in GTEx_v10 folder in T-drive
###############################################################################################################


# Now get normalized count 

dge <- DGEList(counts = cts)

# calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# get normalized counts
normalized_counts <- cpm(dge)
head(normalized_counts[1:10,1:10])
dim(normalized_counts)

# Set function
inv_norm_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# apply inverse normal transform to each gene
final_normalized <- t(apply(normalized_counts, 1, inv_norm_transform))
final_normalized <- as.data.frame(final_normalized)

# add gene annotation back (NO hard-coded 349/348)
final_normalized$GENEID <- rownames(final_normalized)

final_normalized$GENESYMBOL <- expected_counts$GENESYMBOL[match(final_normalized$GENEID,
                                                                expected_counts$GENEID)]
# order columns: gene info first, then samples
final_normalized <- final_normalized[, c("GENESYMBOL", "GENEID", colnames(cts))]
final_normalized <- final_normalized[,c(-2,-3)]
# rename to fit the sctipt 
qnorm.DeseqNC.set <- final_normalized


# ###############################################################################################################
# ###################################  For merge TPM data  ############################
# ###############################################################################################################
# cts <- as.data.frame(cts)
# # add gene annotation back (NO hard-coded 349/348)
# cts$GENEID <- rownames(cts)
# cts$GENESYMBOL <- expected_counts$GENESYMBOL[match(cts$GENEID,
#                                                    expected_counts$GENEID)]
# # order columns: gene info first, then samples
# cts <- cts[, c("GENESYMBOL", "GENEID", colnames(cts)[1:114])]
# qnorm.DeseqNC.set <- cts
# ###############################################################################################################




###############################################################################################################
###############################################################################################################

### Load VcfR 
library(vcfR)
library(biomaRt)

# ==============================================================================
# 1) LOAD VCF AND EXTRACT DATA
# ==============================================================================

#vcf_path <- "/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_500k_uniq.vcf.gz"
#vcf_path <- "/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3.unique.vcf.gz"

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

# As the GT data from GDCconfluence may be flipped and differ from dbSNP, need to adjust

# ==============================================================================
# dbSNP155 (GRCh38) rsID + hg38 strand check + optional genotype flip (OFFLINE)
# Replace your old BIOMART block with this entire section
# ==============================================================================

library(dplyr)
library(GenomicRanges)
library(SNPlocs.Hsapiens.dbSNP155.GRCh38)
library(BSgenome.Hsapiens.UCSC.hg38)

# SNPlocs.Hsapiens.dbSNP155.GRCh38 uses NCBI style: "1", "2", ... not "chr1"
geno_with_pos$CHR <- gsub("^chr", "", geno_with_pos$CHR, ignore.case = TRUE)

# ---- 2) Add rsIDs from dbSNP155 by CHR+POS ----
gr <- GRanges(
  seqnames = geno_with_pos$CHR,
  ranges   = IRanges(start = geno_with_pos$POS, end = geno_with_pos$POS)
)

locs <- snpsByOverlaps(SNPlocs.Hsapiens.dbSNP155.GRCh38, gr)


snp_results <- data.frame(
  CHR       = as.character(seqnames(locs)),
  POS       = start(locs),
  refsnp_id = paste0(as.character(mcols(locs)$RefSNP_id)),
  stringsAsFactors = FALSE
)

# keep one rsID per CHR:POS (dbSNP can return multiple)
snp_results <- snp_results[!duplicated(paste(snp_results$CHR, snp_results$POS, sep=":")), ]


# ---- 3) Merge rsIDs onto your VCF-derived table ----
final_data <- dplyr::left_join(
  geno_with_pos,
  snp_results,
  by = c("CHR", "POS"))

final_data$Clean_SNP_ID <- ifelse(!is.na(final_data$refsnp_id),
                                  final_data$refsnp_id,
                                  final_data$SNP_ID)


# ---- 4) Strand check vs hg38 reference (authoritative) ----
complement_seq <- function(x) chartr("ATCG", "TAGC", x)

final_data$CHR <- paste0('chr',"",final_data$CHR)

hg38_ref <- getSeq(
  BSgenome.Hsapiens.UCSC.hg38,
  names = final_data$CHR,
  start = final_data$POS,
  end   = final_data$POS
)

final_data$hg38_ref_base <- as.character(hg38_ref)

# Only meaningful for simple SNPs (single base REF)
is_simple_snp <- !is.na(final_data$REF) & nchar(final_data$REF) == 1 & final_data$REF %in% c("A","C","G","T")

final_data$Strand_Status <- "Unknown"
final_data$Strand_Status[is_simple_snp & final_data$REF == final_data$hg38_ref_base] <- "Forward"
final_data$Strand_Status[is_simple_snp & complement_seq(final_data$REF) == final_data$hg38_ref_base] <- "Reverse/Flipped"
final_data$Strand_Status[is_simple_snp &
                           !(final_data$REF == final_data$hg38_ref_base) &
                           !(complement_seq(final_data$REF) == final_data$hg38_ref_base)] <- "Mismatch"

# ---- 5) Optional: flip genotype allele strings when Reverse/Flipped (SNP-only) ----
sample_names <- colnames(vcf_genotypes.tmp)

is_snp_gt <- function(gt) !is.na(gt) && grepl("^[ACGT]/[ACGT]$", gt)

for (sample in sample_names) {
  final_data[[sample]] <- ifelse(
    final_data$Strand_Status == "Reverse/Flipped" &
      vapply(final_data[[sample]], is_snp_gt, logical(1)),
    complement_seq(final_data[[sample]]),
    final_data[[sample]]
  )
}

# ---- 6) Final columns (metadata first, then samples) ----
info_cols <- c("CHR", "POS", "SNP_ID", "Clean_SNP_ID", "REF", "ALT",
               "hg38_ref_base", "Strand_Status")

final_data <- final_data[, c(info_cols, sample_names)]



# modified the sample name 
# final_data <- geno_with_pos
# Keep the final data that both REF and ALT are not NA
final_data <- final_data %>%
  dplyr::filter(!is.na(REF) & !is.na(ALT))

colnames(qnorm.DeseqNC.set)
# Apply to specific columns (7 to end)
colnames(final_data)[9:ncol(final_data)] <- sub("^(.+)_\\1$", "\\1", colnames(final_data)[9:ncol(final_data)])



# load sample manifest 
gt_sample_man <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Sample_Genotyping/SR0325-017_1_QC_Report_48202503941.xlsx')
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

FGFR3_qn <- qnorm.DeseqNC.set %>% dplyr::filter(GENESYMBOL == "FGFR3")
FGFR3_qn <- t(FGFR3_qn) %>% as.data.frame()
FGFR3_qn$Sample_Name <- rownames(FGFR3_qn)
FGFR3_qn <- FGFR3_qn[-1,]
FGFR3_qn <- FGFR3_qn[,c(9,1:8)]


final_data_sub <- final_data
# Test it use rsid only 
#final_data_sub <- final_data[grep("rs",final_data$SNP_ID),]
final_data_sub <- final_data_sub[,c(4,9:ncol(final_data_sub))]
final_data_sub <- t(final_data_sub) %>% as.data.frame()
names(final_data_sub) <- final_data_sub[1,]
final_data_sub <- final_data_sub[-1,]

# remove duplicateion 
final_data_sub <- final_data_sub[,!duplicated(names(final_data_sub))]

base_id <- sub("\\.[0-9]+$", "", rownames(final_data_sub))
keep <- !duplicated(base_id)
final_data_sub <- final_data_sub[keep, ]
rownames(final_data_sub) <- base_id[keep]


final_data_sub <- final_data_sub %>% mutate(Sample_Name=rownames(final_data_sub),.before = 1)

FGFR3_Merged_data <- left_join(FGFR3_qn,final_data_sub,by="Sample_Name")

# Load manifest file of sample age sex etc 
gt_sample_man_sub <- gt_sample_man %>% dplyr::select(Sample_Name,Predicted_Sex,AFR,EUR,ASN,Ancestry,`Sample Pass QC`)

# Load RIN score file
RIN <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Sample_Genotyping/DNA_RNA_Log_for_Tissues.xlsx',sheet = 7)
RIN$Sample_Name <- paste0(RIN$`Sample Type...2`,"_",RIN$Sample...1)
RIN <- RIN[,c(36,12)]
names(RIN) <- c("Sample_Name","RIN_score")
# replace NA RIN score to the lowest 
# Method 1: Simple replacement
RIN$RIN_score[is.na(RIN$RIN_score)] <- min(RIN$RIN_score, na.rm = TRUE)


gt_sample_man_sub <- gt_sample_man_sub %>%
  left_join(RIN, by = "Sample_Name") %>%
  relocate(RIN_score, .after = Sample_Name)

gt_sample_man_sub <- gt_sample_man_sub[!duplicated(gt_sample_man_sub$Sample_Name), ]

final_FGFR3 <- left_join(gt_sample_man_sub,FGFR3_Merged_data,by = 'Sample_Name')

# Done save 
# write.table(final_FGFR3,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Parse_scRNA-seq/Sample_Genotyping/FGFR3_isoform_TMM_INT.snp.csv',row.names = F,quote = F,sep = ',')
# write.table(df_fgfr,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_isoform_TMM_INT.snp_500k.csv',row.names = F,quote = F,sep = ',')

###########################################################
## lm analysis ref from BL project #######
#################################################################################################################################################################################
detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
detach("package:ensembldb", unload = TRUE)

library(dplyr)
library(qvalue)
library(lmPerm)

df_fgfr <- final_FGFR3

# names(df_fgfr) <- sub("(:[^:]+){2}$", "", names(df_fgfr))
# names(df_fgfr)[grep("^4",names(df_fgfr))] <-  paste0('chr',"",names(df_fgfr)[grep("^4",names(df_fgfr))])
# names(df_fgfr)[grep("^chr4",names(df_fgfr))] <- gsub(":","_",names(df_fgfr)[grep("^chr4",names(df_fgfr))])

df_fgfr <- df_fgfr %>% 
  filter(!is.na(Predicted_Sex) & Predicted_Sex != "") %>%
  filter(!is.na(Ancestry) & Ancestry != "") 


#### check correct asign factor or chr or numeric

df_fgfr$Predicted_Sex <- as.factor(df_fgfr$Predicted_Sex)
df_fgfr$Ancestry <-as.factor(df_fgfr$Ancestry)


# ---- Make simple SNP -> (hg38 REF, hg38-aligned ALT) lookup ----
complement_seq <- function(x) chartr("ATCG", "TAGC", x)

allele_lut <- final_data[, c("Clean_SNP_ID", "hg38_ref_base", "ALT", "Strand_Status")]
allele_lut <- allele_lut[!duplicated(allele_lut$Clean_SNP_ID), ]

allele_lut$ALT_hg38 <- allele_lut$ALT
flip <- allele_lut$Strand_Status == "Reverse/Flipped" &
  !is.na(allele_lut$ALT) & nchar(allele_lut$ALT) == 1
allele_lut$ALT_hg38[flip] <- complement_seq(allele_lut$ALT[flip])


#### Make SNP dosage as 0/1/2 = number of ALT alleles (hg38 aligned)

for (i in names(df_fgfr)[17:ncol(df_fgfr)]) {
  
  # Fallback map for VCF REF/ALT (used only if Strand_Status == "Mismatch")
  vcf_map <- final_data[, c("Clean_SNP_ID", "REF", "ALT")]
  vcf_map <- vcf_map[!duplicated(vcf_map$Clean_SNP_ID), ]
  
  # i <- "rs1665366"
  n1 <- paste0(i, "_add")
  
  idx <- match(i, allele_lut$Clean_SNP_ID)
  if (is.na(idx)) next
  
  REF_use <- allele_lut$hg38_ref_base[idx]
  ALT_use <- allele_lut$ALT_hg38[idx]
  
  # If mismatch, use original VCF REF/ALT instead
  if (allele_lut$Strand_Status[idx] == "Mismatch") {
    idx2 <- match(i, vcf_map$Clean_SNP_ID)
    if (!is.na(idx2)) {
      REF_use <- vcf_map$REF[idx2]
      ALT_use <- vcf_map$ALT[idx2]
    }
  }
  
  if (is.na(REF_use) || is.na(ALT_use) || REF_use == ALT_use) next
  
  df_fgfr[[n1]] <- dplyr::case_when(
    df_fgfr[[i]] == paste0(REF_use, "/", REF_use) ~ 0,
    df_fgfr[[i]] %in% c(paste0(REF_use, "/", ALT_use), paste0(ALT_use, "/", REF_use)) ~ 1,
    df_fgfr[[i]] == paste0(ALT_use, "/", ALT_use) ~ 2,
    TRUE ~ NA_real_
  )
}



df_fgfr <- df_fgfr[,c(names(df_fgfr)[1:16],setdiff(names(df_fgfr),names(df_fgfr)[1:16]) %>% sort())]


#############################

# Can make metadata list 

megdata <- list()

megdata[["FGFR3"]] <- df_fgfr

inputdf <- megdata$FGFR3  


####################################################################
####################################################################
################# NOW WE DO lm() ###################################
####################################################################
####################################################################

# use the saved file 
#inputdf <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_isoform_TMM_INT.snp_500k.csv')

inputdf <- df_fgfr

str(inputdf)

inputdf$Predicted_Sex <- factor(inputdf$Predicted_Sex)
inputdf$Ancestry <- factor(inputdf$Ancestry)

df.out <- data.frame()
df_count.summary <- data.frame()
# 



for (i in grep("_add",names(inputdf),value = T)){
  # i <- "rs62286992_add"
  
  snp_id <- gsub("_add", "", i)   # "rs62286992"
  
  m <- allele_lut[match(snp_id, allele_lut$Clean_SNP_ID), ]
  REF <- m$hg38_ref_base
  ALT <- m$ALT_hg38
  Strand_Status <- m$Strand_Status
  
  
  
  
  function.names <- c("max", "min", "mean", "median", "sd")
  
  for ( k in function.names ) {
    # k <- function.names[3]
    # the all the TPM in all tissue
    # consider dealing with NA
    
    # & inputdf[[i]] != 0
    # names(inputdf[6:135]
    
    # remove SNPs that are NA and empty
    df.tmp <- inputdf[,c(names(inputdf[4]), gsub("_add", "", i))] %>%
      filter(!is.na(inputdf[[gsub("_add", "", i)]]) & inputdf[[gsub("_add", "", i)]] != ""  ) %>%
      group_by(get(noquote(gsub("_add", "", i)))) %>%
      dplyr::summarise(across(where(is.numeric), get(k), na.rm=TRUE)) %>%
      data.frame()
    
    df.tmp$variable <- i
    df.tmp$stat <- k
    
    df.tmp <- df.tmp %>% dplyr::select(c("stat", "variable", everything())) 
    names(df.tmp)[3] <- "genotype"
    
    df_count.summary <- rbind(df_count.summary, df.tmp)
    df_count.summary[] <- lapply(df_count.summary, function(x) if(is.numeric(x)) round(x, 2) else x)
    
    
    
    
  }
  
  
  
  # names(inputdf[6:135]
  
  for(j in names(inputdf)[c(10,16,11)]){
    
    # j <- "ENST00000481110"
    
    
    # Filter the data to exclude rows where the SNP dosage is 0
    filtered_0_1_2_only <- inputdf %>% filter(inputdf[[i]] %in% c(0, 1, 2))
    
    # counting GT of how many sample 
    GT.count <- filtered_0_1_2_only
    
    GT.count.0a <- GT.count[ which(GT.count[,i] == 0), ]
    GT.count.0b <- as.character(unique(GT.count.0a[gsub("_add", "", i)]))
    GT.count.1a <- GT.count[ which(GT.count[,i] == 1), ]
    GT.count.1b <- as.character(unique(GT.count.1a[gsub("_add", "", i)]))
    GT.count.2a <- GT.count[ which(GT.count[,i] == 2), ]
    GT.count.2b <- as.character(unique(GT.count.2a[gsub("_add", "", i)]))
    
    # dynamically generate formula
    fmla_un <- as.formula(paste0(j, "~" , i))
    fmla_adj1 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex +  Ancestry + RIN_score"))
    fmla_adj2 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex + AFR + EUR + ASN + RIN_score"))
    fmla_adj_int1 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex + RIN_score + Ancestry +",i,paste0("*Predicted_Sex"))) 
    fmla_adj_int2 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex + RIN_score + AFR + EUR + ASN + ",i,paste0("*Predicted_Sex")))
    
    
    
    ######################################################
    ######### make 0 as 1, not sure we need here #########
    # filtered_0_1_2_only <- filtered_0_1_2_only %>% mutate(!!i := ifelse(.[[i]] == 0, 1, .[[i]]))
    ######################################################
    
    # fit lm model
    fit_un <- lm(fmla_un, filtered_0_1_2_only)
    fit_adj1 <- lm(fmla_adj1, filtered_0_1_2_only)
    fit_adj2 <- lm(fmla_adj2, filtered_0_1_2_only)
    fit_adj_int1 <- lm(fmla_adj_int1, filtered_0_1_2_only)
    fit_adj_int2 <- lm(fmla_adj_int2, filtered_0_1_2_only)
    
    ########################################
    # run permutation #
    ########################################
    
    #prm <- lmp(fmla_un, data = filtered_0_1_2_only, perm = "Prob")
    #prm_adj <- lmp(fmla_adj, data = filtered_0_1_2_only, perm = "Prob")
    # Calculate the number of non-zero samples in j where i is not NA
    # use this since the lm() is done based on the sample inculde GT
    
    non_zero_count <- sum(filtered_0_1_2_only[[j]] != 0 & !is.na(filtered_0_1_2_only[[j]]))
    
    # non_zero_count <- sum(inputdf[[j]] != 0 & !is.na(inputdf[[j]]))
    
    # Calculate mean values for SNP dosages 1 and 2
    # mean_value_1 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 1) %>% pull(j), na.rm = TRUE),2) 
    # mean_value_2 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 2) %>% pull(j), na.rm = TRUE),2) 
    
    # create temporary data frame
    df.lm <- data.frame(
      
      snp = snp_id,
      Strand_Status = Strand_Status,
      REF = REF,
      ALT = ALT,
      variable = j,
      sample_used_in_analysis = length(filtered_0_1_2_only[[j]]),
      sample_over_zero = non_zero_count,
      perc_sample_over_zero = round((non_zero_count / length(filtered_0_1_2_only[[j]])) * 100,2),
      
      geno_0 = GT.count.0b,
      geno_1 = GT.count.1b,
      geno_2 = GT.count.2b,
      
      n_0 = nrow(GT.count.0a),
      n_1 = nrow(GT.count.1a),
      n_2 = nrow(GT.count.2a),
      
      #mean_value_0_and_1=mean_value_1,
      #mean_value_2=mean_value_2,
      
      model_un = paste0(formula(fit_un)[2],"~",formula(fit_un)[3]),
      model_adj1 = paste0(formula(fit_adj1)[2],"~",formula(fit_adj1)[3]),
      model_adj2 = paste0(formula(fit_adj2)[2],"~",formula(fit_adj2)[3]),
      model_adj_int1 = paste0(formula(fit_adj_int1)[2],"~",formula(fit_adj_int1)[3]),
      model_adj_int2 = paste0(formula(fit_adj_int2)[2],"~",formula(fit_adj_int2)[3]),
      
      # add lm() result
      
      #p.prm = tryCatch({coef(summary(prm))[2,3]}, error=function(e){
      #  print(paste0('error of NA'))
      #  return(NA)}),
      
      # Main effects (wrapped in tryCatch to prevent loop breaks)
      p.value_un = tryCatch(coef(summary(fit_un))[2,4]%>% round(.,3), error = function(e) NA),
      beta_un = tryCatch(round(coef(summary(fit_un))[2,1], 3), error = function(e) NA),
      
      #StdEr = round(coef(summary(fit_un))[2,2], digits = 4),
      
      p.value_adj1 = coef(summary(fit_adj1))[2,4]%>% round(.,3),
      beta_adj1 = coef(summary(fit_adj1))[2,1]%>% round(.,3),
      
      # StdEr_adj_sex_age = round(coef(summary(fit_adj))[2,2], digits = 4),
      # p.prm_adj_sex_age = coef(summary(prm_adj))[2,3],
      
      p.value_adj2 = coef(summary(fit_adj2))[2,4] %>% round(.,3),
      beta_adj2 = coef(summary(fit_adj2))[2,1]%>% round(.,3),
      
      p.value_adj_int1 = tryCatch(coef(summary(fit_adj_int1))[grep(paste0("^",i, ":"), rownames(coef(summary(fit_adj_int1)))), "Pr(>|t|)"]%>% round(.,3),error = function(e) NA)[1],
      p.value_adj_int2 = tryCatch(coef(summary(fit_adj_int2))[grep(paste0("^",i, ":"), rownames(coef(summary(fit_adj_int2)))), "Pr(>|t|)"]%>% round(.,3),error = function(e) NA)[1]
      
    )
    
    # bind rows of temporary data frame to the results data frame
    
    df.out <- rbind(df.out, df.lm)
    
    
  }
}

