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

#####################################################################################
###################################  For merge TPM data  ############################
#####################################################################################
cts <- as.data.frame(cts)
# add gene annotation back (NO hard-coded 349/348)
cts$GENEID <- rownames(cts)
cts$GENESYMBOL <- expected_counts$GENESYMBOL[match(cts$GENEID,
                                                   expected_counts$GENEID)]
# order columns: gene info first, then samples
cts <- cts[, c("GENESYMBOL", "GENEID", colnames(cts)[1:114])]
qnorm.DeseqNC.set <- cts

###############################################################################################################
###############################################################################################################

### Load VcfR 
library(vcfR)
library(biomaRt)

# ==============================================================================
# 1) LOAD VCF AND EXTRACT DATA
# ==============================================================================

vcf_path <- "/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_500k_uniq.vcf.gz"
#vcf_path <- "/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3.500k.SUBs.sub.V2.vcf.gz"

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

###****NEW

is_simple_snp <- !is.na(final_data$REF) & nchar(final_data$REF) == 1 & final_data$REF %in% c("A","C","G","T") &
  !is.na(final_data$ALT) & nchar(final_data$ALT) == 1 & final_data$ALT %in% c("A","C","G","T")

ref_match <- final_data$REF == final_data$hg38_ref_base
alt_match <- final_data$ALT == final_data$hg38_ref_base
ref_comp_match <- complement_seq(final_data$REF) == final_data$hg38_ref_base
alt_comp_match <- complement_seq(final_data$ALT) == final_data$hg38_ref_base

final_data$Strand_Status <- "Unknown"
final_data$Strand_Status[is_simple_snp & (ref_match | alt_match)] <- "Forward"
final_data$Strand_Status[is_simple_snp & (ref_comp_match | alt_comp_match)] <- "Reverse/Flipped"
final_data$Strand_Status[is_simple_snp & !(ref_match | alt_match | ref_comp_match | alt_comp_match)] <- "Mismatch"

####***NEW

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
if 
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
RIN$Sample_Name[23] <- "VR_35_D"


gt_sample_man_sub <- gt_sample_man_sub %>%
  left_join(RIN, by = "Sample_Name") %>%
  relocate(RIN_score, .after = Sample_Name)

gt_sample_man_sub <- gt_sample_man_sub[!duplicated(gt_sample_man_sub$Sample_Name), ]

final_FGFR3 <- left_join(gt_sample_man_sub,FGFR3_Merged_data,by = 'Sample_Name')

###########################################################
## lm analysis ref from BL project #######
#################################################################################################################################################################################
detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
detach("package:ensembldb", unload = TRUE)

library(dplyr)
library(qvalue)
library(lmPerm)

df_fgfr <- final_FGFR3

df_fgfr <- df_fgfr %>% 
  filter(!is.na(Predicted_Sex) & Predicted_Sex != "") %>%
  filter(!is.na(Ancestry) & Ancestry != "") 


#### check correct asign factor or chr or numeric

df_fgfr$Predicted_Sex <- as.factor(df_fgfr$Predicted_Sex)
df_fgfr$Ancestry <-as.factor(df_fgfr$Ancestry)


# ---- Make simple SNP -> (hg38 REF, hg38-aligned ALT) lookup ----
complement_seq <- function(x) chartr("ACGT", "TGCA", x)

allele_lut <- final_data[, c("Clean_SNP_ID", "hg38_ref_base", "REF", "ALT", "Strand_Status")]
allele_lut <- allele_lut[!duplicated(allele_lut$Clean_SNP_ID), ]

ref2 <- allele_lut$REF
alt2 <- allele_lut$ALT

flip <- allele_lut$Strand_Status == "Reverse/Flipped" &
  !is.na(ref2) & nchar(ref2) == 1 &
  !is.na(alt2) & nchar(alt2) == 1

ref2[flip] <- complement_seq(ref2[flip])
alt2[flip] <- complement_seq(alt2[flip])

allele_lut$ALT_hg38 <- NA_character_

m_ref <- ref2 == allele_lut$hg38_ref_base
m_alt <- alt2 == allele_lut$hg38_ref_base

allele_lut$ALT_hg38[m_ref] <- alt2[m_ref]
allele_lut$ALT_hg38[m_alt] <- ref2[m_alt]


# Fallback map for VCF REF/ALT (used only if Strand_Status == "Mismatch")
vcf_map <- final_data[, c("Clean_SNP_ID", "REF", "ALT")]
vcf_map <- vcf_map[!duplicated(vcf_map$Clean_SNP_ID), ]

########*NEW Now load LD_proxy result so that we can chkec risk protect allele
########*
## Get the proxy from LD-link R package, find the SNP that exist in the proxy data
# Use LDlinkR to extract SNP proxy
library(LDlinkR)
library(stringr)

# EUR hg38 coord 
my_proxies <- LDproxy(snp = "rs2896518", 
                      pop = "EUR", 
                      r2d = "r2", 
                      token = "c22a071b5bda",genome_build = "grch38_high_coverage",win_size = 500000
)

my_proxies$Coord <- gsub(':','_',my_proxies$Coord)

my_proxies <- my_proxies %>%
  mutate(
    risk.Correlated    = str_extract(Correlated_Alleles, "(?<=A=)[^,]+"),
    protect.Correlated = str_extract(Correlated_Alleles, "(?<=G=).+")
  )  %>% relocate(risk.Correlated,protect.Correlated,.before = FORGEdb)

# Select the SNPs that exist in the 
Prox_SNP_list <- my_proxies$RS_Number %>% grep('rs',.,value = T)

df_fgfr_prox <- df_fgfr[, c(names(df_fgfr)[1:16],
                            intersect(names(df_fgfr), Prox_SNP_list))]

# Now set rsik as 2 and protect as 0
#### Make SNP dosage as 0/1/2 = number of ALT alleles (hg38 aligned)

# ---- Create risk allele lookup ----
risk_lut <- my_proxies %>%
  dplyr::select(RS_Number, risk.Correlated, protect.Correlated) %>%
  dplyr::distinct(RS_Number, .keep_all = TRUE)

snp_cols <- names(df_fgfr_prox)[17:ncol(df_fgfr_prox)]

for (i in snp_cols) {
  #i <- "rs798727"
  #i <- "rs2402494"
  #i <- "rs67275758"
  n1 <- paste0(i, "_add")
  
  ridx <- match(i, risk_lut$RS_Number)
  RISK_use <- if (!is.na(ridx)) risk_lut$risk.Correlated[ridx] else NA_character_
  PROT_use <- if (!is.na(ridx)) risk_lut$protect.Correlated[ridx] else NA_character_
  
  # Fallback alleles from your hg38-aligned lookup
  aidx <- match(i, allele_lut$Clean_SNP_ID)
  if (is.na(aidx)) next
  
  REF_use <- allele_lut$hg38_ref_base[aidx]
  ALT_use <- allele_lut$ALT_hg38[aidx]
  
  if (allele_lut$Strand_Status[aidx] == "Mismatch") {
    vidx <- match(i, vcf_map$Clean_SNP_ID)
    if (!is.na(vidx)) {
      REF_use <- vcf_map$REF[vidx]
      ALT_use <- vcf_map$ALT[vidx]
    }
  }
  
  if (is.na(REF_use) || is.na(ALT_use)) next
  
  gt <- df_fgfr_prox[[i]]
  
  # --- Case 1: we have risk/protect from LDproxy -> risk dosage ---
  if (!is.na(RISK_use) && !is.na(PROT_use)) {
    df_fgfr_prox[[n1]] <- dplyr::case_when(
      gt == paste0(PROT_use, "/", PROT_use) ~ 0,
      gt %in% c(paste0(RISK_use, "/", PROT_use), paste0(PROT_use, "/", RISK_use)) ~ 1,
      gt == paste0(RISK_use, "/", RISK_use) ~ 2,
      TRUE ~ NA_real_
    )
    next
  }
  
  # --- Case 2: missing risk/protect -> ALT dosage fallback ---
  df_fgfr_prox[[n1]] <- dplyr::case_when(
    gt == paste0(REF_use, "/", REF_use) ~ 0,
    gt %in% c(paste0(REF_use, "/", ALT_use), paste0(ALT_use, "/", REF_use)) ~ 1,
    gt == paste0(ALT_use, "/", ALT_use) ~ 2,
    TRUE ~ NA_real_
  )
}

df_fgfr_prox <- df_fgfr_prox[,c(names(df_fgfr_prox)[1:16],setdiff(names(df_fgfr_prox),names(df_fgfr_prox)[1:16]) %>% sort())]

########*NEW
########*
########*


####################################################################
####################################################################
################# NOW WE DO lm() ###################################
####################################################################
####################################################################

# use the saved file 
#inputdf <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_isoform_TMM_INT.snp_500k_v3.csv')
# inputdf <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/FGFR3_isoform_TPM.snp_500k_v3.csv')



data_list <- list()
data_list$Vict <- df_fgfr_prox[grep("VR",df_fgfr_prox$Sample_Name),]
data_list$Tissue_piece <- df_fgfr_prox[grep("DS",df_fgfr_prox$Sample_Name),]


# inputdf <- data_list$Vict

names(inputdf)[c(10,16,11)] <- c("FGFR3_IIIb","FGFR3_IIIc","FGFR3_IIIs")

str(inputdf)


inputdf$Predicted_Sex <- factor(inputdf$Predicted_Sex)
inputdf$Ancestry <- factor(inputdf$Ancestry)
inputdf[, 9:16] <- lapply(inputdf[, 9:16], as.numeric)

df.out <- data.frame()
df_count.summary <- data.frame()

for (i in grep("_add",names(inputdf),value = T)){
  # i <- "rs4865455_add"
  # 
  snp_id <- gsub("_add", "", i)   # "rs62286992"
  
  m <- allele_lut[match(snp_id, allele_lut$Clean_SNP_ID), ]
  
  # If mismatch, use original VCF REF/ALT instead
  if (m$Strand_Status == "Mismatch") {
    
    idx2 <- match(snp_id, vcf_map$Clean_SNP_ID)
    REF <- vcf_map$REF[idx2]
    ALT <- vcf_map$ALT[idx2]} else {
      
      REF <- m$hg38_ref_base
      ALT <- m$ALT_hg38
      
    }
  
  
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
    
    # j <- "FGFR3-IIIb"
    
    
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
    
    #non_zero_count <- sum(filtered_0_1_2_only[[j]] != 0 & !is.na(filtered_0_1_2_only[[j]]))
    # non_zero_count <- sum(inputdf[[j]] != 0 & !is.na(inputdf[[j]]))
    
    
    
    # create temporary data frame
    df.lm <- data.frame(
      
      snp = snp_id,
      Strand_Status = Strand_Status,
      REF = REF,
      ALT = ALT,
      variable = j,
      
      #sample_used_in_analysis = length(filtered_0_1_2_only[[j]]),
      #sample_over_zero = non_zero_count,
      # perc_sample_over_zero = round((non_zero_count / length(filtered_0_1_2_only[[j]])) * 100,2),
      
      geno_0 = GT.count.0b,
      geno_1 = GT.count.1b,
      geno_2 = GT.count.2b,
      
      n_0 = nrow(GT.count.0a),
      n_1 = nrow(GT.count.1a),
      n_2 = nrow(GT.count.2a),
      
      # Calculate mean values for SNP dosages 1 and 2
      
      #TPM_mean = round(mean(filtered_0_1_2_only$FGFR3_total_TPM),2),
      exp_geno_0 = round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 0) %>% pull(j), na.rm = TRUE),2),
      exp_geno_1 = round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 1) %>% pull(j), na.rm = TRUE),2),
      exp_geno_2 = round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 2) %>% pull(j), na.rm = TRUE),2),
      
      #mean_value_0_and_1=mean_value_1,
      #mean_value_2=mean_value_2,
      
      #model_un = paste0(formula(fit_un)[2],"~",formula(fit_un)[3]),
      #model_adj1 = paste0(formula(fit_adj1)[2],"~",formula(fit_adj1)[3]),
      #model_adj2 = paste0(formula(fit_adj2)[2],"~",formula(fit_adj2)[3]),
      #model_adj_int1 = paste0(formula(fit_adj_int1)[2],"~",formula(fit_adj_int1)[3]),
      #model_adj_int2 = paste0(formula(fit_adj_int2)[2],"~",formula(fit_adj_int2)[3]),
      
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

#######

# read in house data
Indf <- inputdf

Our_target <- names(Indf) %>% grep('chr4_[0-9+]+$|rs[0-9]+$',.,value = T) %>% as.data.frame()
names(Our_target) <- "Coord"

target_for_rsid <- Our_target %>% filter(grepl("^rs[0-9]+$", Coord))
names(target_for_rsid) <- "RS_Number"
cb_for_rsid <- inner_join(target_for_rsid,my_proxies,by = 'RS_Number')

# Extract the lm result to extract the SNPs that is in proxy data

# res_lm <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/lm_res/')

res_lm <- df.out

names(res_lm)[1] <- "RS_Number"

res_lm_rsONLY <- res_lm[(grepl('rs',res_lm$RS_Number)),]

COMBINE_rs <- inner_join(cb_for_rsid,res_lm_rsONLY,by = 'RS_Number')

library(stringr)

COMBINE_rs <- COMBINE_rs %>%
  mutate(
    risk.Correlated    = str_extract(Correlated_Alleles, "(?<=A=)[^,]+"),
    protect.Correlated = str_extract(Correlated_Alleles, "(?<=G=).+")
  )  %>% relocate(risk.Correlated,protect.Correlated,.before = variable)




