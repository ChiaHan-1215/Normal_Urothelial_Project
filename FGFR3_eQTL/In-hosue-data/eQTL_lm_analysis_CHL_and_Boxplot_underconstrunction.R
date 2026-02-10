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
#write.table('')
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

# # ==============================================================================
# # 2) BIOMART QUERY (hg38)
# # ==============================================================================
# query_chr <- gsub("chr", "", geno_with_pos$CHR, ignore.case = TRUE)
# coords <- paste(query_chr, geno_with_pos$POS, geno_with_pos$POS, sep = ":")
# 
# snpMart <- useEnsembl(biomart = "snp", dataset = "hsapiens_snp")
# 
# snp_results_raw <- getBM(
#   attributes = c('refsnp_id', 'chrom_start', 'allele'),
#   filters = 'chromosomal_region', 
#   values = coords, 
#   mart = snpMart
# )
# 
# # Filter for exact coordinate matches
# snp_results <- snp_results_raw[snp_results_raw$chrom_start %in% geno_with_pos$POS, ]
# 
# # ==============================================================================
# # 3) MERGE AND STRAND DETECTION
# # ==============================================================================
# final_data <- merge(
#   geno_with_pos, 
#   snp_results, 
#   by.x = "POS", 
#   by.y = "chrom_start",
#   all.x = TRUE
# )
# 
# # Function to get the complement (A <-> T, C <-> G)
# complement_seq <- function(sequence) {
#   # Handles single bases or strings like "A/C"
#   chartr("ATCG", "TAGC", sequence)
# }
# 
# # Identify Strand Status
# final_data$Strand_Status <- sapply(1:nrow(final_data), function(i) {
#   vcf_ref <- final_data$REF[i]
#   ens_alleles <- final_data$allele[i]
#   if (is.na(ens_alleles)) return("Unknown")
#   
#   if (grepl(vcf_ref, ens_alleles)) {
#     return("Forward")
#   } else if (grepl(complement_seq(vcf_ref), ens_alleles)) {
#     return("Reverse/Flipped")
#   } else {
#     return("Mismatch")
#   }
# })
# 
# # ==============================================================================
# # 4) AUTOMATIC GENOTYPE UNF-FLIPPING (The Fix)
# # ==============================================================================
# # This section iterates through all sample columns and complements the alleles
# # ONLY if the Strand_Status is "Reverse/Flipped".
# 
# sample_names <- colnames(vcf_genotypes.tmp)
# 
# for (sample in sample_names) {
#   # Apply complement only where flipped
#   # sample <- sample_names[2]
#   final_data[[sample]] <- ifelse(
#     final_data$Strand_Status == "Reverse/Flipped",
#     complement_seq(final_data[[sample]]), 
#     final_data[[sample]]
#   )
# }
# 
# # ==============================================================================
# # 5) FINAL CLEANUP AND EXPORT
# # ==============================================================================
# # Update IDs to official rsIDs
# final_data$Clean_SNP_ID <- ifelse(
#   !is.na(final_data$refsnp_id), final_data$refsnp_id, final_data$SNP_ID
# )
# 
# # Rename Ensembl column
# colnames(final_data)[colnames(final_data) == "allele"] <- "Ensembl_Alleles_Plus_Strand"
# 
# # Organize columns: Metadata first, then Samples
# info_cols <- c("CHR", "POS","SNP_ID","Clean_SNP_ID", "REF", "ALT", 
#                "Ensembl_Alleles_Plus_Strand", "Strand_Status")
# final_data <- final_data[, c(info_cols, sample_names)]
# 
# 
# 

# modified the sample name 
final_data <- geno_with_pos
# Keep the final data that both REF and ALT are not NA
final_data <- final_data %>%
  dplyr::filter(!is.na(REF) & !is.na(ALT))

colnames(qnorm.DeseqNC.set)
# Apply to specific columns (7 to end)
colnames(final_data)[6:ncol(final_data)] <- sub("^(.+)_\\1$", "\\1", colnames(final_data)[6:ncol(final_data)])



# load sample manifest 
gt_sample_man <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Sample_Genotyping/SR0325-017_1_QC_Report_48202503941.xlsx')
gt_sample_man$Sample_Name <- gsub(" ","_",gt_sample_man$Sample_Name)
gt_sample_man$Sample_Name <- gsub("_DNA","",gt_sample_man$Sample_Name)

# cahnge name 

guide <- setNames(gt_sample_man$Sample_Name,gt_sample_man$Sample_ID)

old_names <- colnames(final_data)[6:ncol(final_data)]

# Replace only where a match exists in guide
new_names <- ifelse(
  old_names %in% names(guide),
  guide[old_names],
  old_names
)

colnames(final_data)[6:ncol(final_data)] <- new_names


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
final_data_sub <- final_data_sub[,c(3,6:125)]
final_data_sub <- t(final_data_sub) %>% as.data.frame()
names(final_data_sub) <- final_data_sub[1,]
final_data_sub <- final_data_sub[-1,]

# remove duplicateion 
final_data_sub <- final_data_sub[,!duplicated(names(final_data_sub))]

final_data_sub <- final_data_sub %>% mutate(Sample_Name=rownames(final_data_sub),.before = 1)

FGFR3_Merged_data <- left_join(FGFR3_qn,final_data_sub,by="Sample_Name")

# Load manifest file of sample age sex etc 
gt_sample_man_sub <- gt_sample_man %>% dplyr::select(Sample_Name,Predicted_Sex,AFR,EUR,ASN,Ancestry,`Sample Pass QC`)

# Load RIN score file
RIN <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Sample_Genotyping/DNA_RNA_Log_for_Tissues.xlsx',sheet = 7)
RIN$Sample_Name <- paste0(RIN$`Sample Type...2`,"_",RIN$Sample...1)
RIN <- RIN[,c(36,12)]
names(RIN) <- c("Sample_Name","RIN_score")

gt_sample_man_sub <- gt_sample_man_sub %>%
  left_join(RIN, by = "Sample_Name") %>%
  relocate(RIN_score, .after = Sample_Name)


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

names(df_fgfr) <- sub("(:[^:]+){2}$", "", names(df_fgfr))
names(df_fgfr)[grep("^4",names(df_fgfr))] <-  paste0('chr',"",names(df_fgfr)[grep("^4",names(df_fgfr))])
names(df_fgfr)[grep("^chr4",names(df_fgfr))] <- gsub(":","_",names(df_fgfr)[grep("^chr4",names(df_fgfr))])


df_fgfr <- df_fgfr %>% 
  filter(!is.na(Predicted_Sex) & Predicted_Sex != "") %>%
  filter(!is.na(Ancestry) & Ancestry != "") 


#### check correct asign factor or chr or numeric

df_fgfr$Predicted_Sex <- as.factor(df_fgfr$Predicted_Sex)
df_fgfr$Ancestry <-as.factor(df_fgfr$Ancestry)


#### Make SNP as 0,1,2 

# Iterate through SNP columns
for (i in names(df_fgfr)[12:79]) {
  
  #i <- "4:1763224:G:C"
  n1 <- paste0(i, "_add")
  
  # 1. Create frequency table and ensure it's a clean data frame
  tp1 <- table(df_fgfr[[i]], useNA = "no")
  if (length(tp1) == 0) next # Skip if column is empty
  
  tp1 <- as.data.frame(tp1, stringsAsFactors = FALSE)
  colnames(tp1) <- c("Genotype", "Freq")
  
  # 2. Split Genotype into Alleles (e.g., "G/G" -> "G", "G")
  tp1 <- tp1 %>%
    separate(Genotype, into = c("A1", "A2"), sep = "/", 
             fill = "right", extra = "merge", remove = FALSE)
  
  # 3. Identify Homozygotes (where Allele 1 == Allele 2)
  # We handle the case where A2 might be NA (e.g. if genotype was just "G")
  SAMEGT <- tp1 %>% 
    filter(!is.na(A1) & !is.na(A2) & A1 == A2) %>%
    arrange(Freq) # Sort by frequency (Ascending: Minor first)
  
  # 4. Identify Heterozygotes (anything not in SAMEGT)
  tp_left <- anti_join(tp1, SAMEGT, by = "Genotype")
  
  # 5. Recode based on the number of homozygotes found
  if (nrow(SAMEGT) == 2) {
    # Typical case: 2 homozygotes (0 and 2) and 1 heterozygote (1)
    recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0; '", 
                         SAMEGT$Genotype[2], "'=2; '", 
                         tp_left$Genotype[1], "'=1")
    df_fgfr[[n1]] <- car::recode(df_fgfr[[i]], recode_str)
    
  } else if (nrow(SAMEGT) == 1) {
    # Case with only 1 homozygote (recode to 0) and potentially 1 heterozygote (recode to 1)
    if (nrow(tp_left) >= 1) {
      recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0; '", tp_left$Genotype[1], "'=1")
    } else {
      recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0")
    }
    df_fgfr[[n1]] <- car::recode(df_fgfr[[i]], recode_str)
  }
}


df_fgfr <- df_fgfr[,c(names(df_fgfr)[1:16],setdiff(names(df_fgfr),names(df_fgfr)[1:16]) %>% sort())]

#############################

# Can make metadata list 

megdata <- list()

megdata[["FGFR3"]] <- df_fgfr
#megdata[["TMB_SBS"]] <- mt_ug


####################################################################
####################################################################
################# NOW WE DO lm() ###################################
####################################################################
####################################################################


# ####################################################################
### TPM/quantile (Z-score) normailzed ####
# ####################################################################


inputdf <- megdata$FGFR3  

# 
df.out <- data.frame()
df_count.summary <- data.frame()
# 




for (i in grep("_add",names(inputdf),value = T)){
  # i <- "4:1768718_add"
  # i <- "rs111457485_add" 
  # i <- "rs16997913_add"
  
  
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
  
  for(j in names(inputdf[4:11])){
    
    # j <- "FGFR3"
    # j <- names(inputdf[8:16])[1]
    # j <- names(inputdf[6:135])[20]
    
    # 3 model, all sample size should be same, add Permutation column
    
    # unadjust: TPM ~ SNP 
    # adjust: TPM ~ SNP + age + sex
    # adjust for EBV: TPM ~ SNP + age + sex +EBV
    
    # for subsets
    # dynamically generate formula
    fmla_un <- as.formula(paste0(j, "~" , i))
    fmla_adj <- as.formula(paste0(j, "~" , i, " + Predicted_Sex + Ancestry "))
    #fmla_adj_ebv <- as.formula(paste0(j, "~" , i, " + Age_final + sex_final + EBV_final"))
    
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
    
    ######################################################
    ######### make 0 as 1, not sure we need here #########
    # filtered_0_1_2_only <- filtered_0_1_2_only %>% mutate(!!i := ifelse(.[[i]] == 0, 1, .[[i]]))
    ######################################################
    
    # fit lm model
    fit_un <- lm(fmla_un, filtered_0_1_2_only)
    fit_adj <- lm(fmla_adj, filtered_0_1_2_only)
    
    
    
    ########################################
    # run permutation #
    ########################################
    
    prm <- lmp(fmla_un, data = filtered_0_1_2_only, perm = "Prob")
    prm_adj <- lmp(fmla_adj, data = filtered_0_1_2_only, perm = "Prob")
    
    # Calculate the number of non-zero samples in j where i is not NA
    
    # use this since the lm() is done based on the sample inculde GT
    non_zero_count <- sum(filtered_0_1_2_only[[j]] != 0 & !is.na(filtered_0_1_2_only[[j]]))
    
    #non_zero_count <- sum(inputdf[[j]] != 0 & !is.na(inputdf[[j]]))
    
    # Calculate mean values for SNP dosages 1 and 2
    mean_value_1 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 1) %>% pull(j), na.rm = TRUE),2) 
    mean_value_2 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 2) %>% pull(j), na.rm = TRUE),2) 
    
    
    
    
    # create temporary data frame
    df.lm <- data.frame(
      
      snp = i,
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
      
      
      
      # add lm() result
      p.value = tryCatch({coef(summary(fit_un))[2,4]}, error=function(e){
        print(paste0('error of NA'))
        return(NA)}) ,
      
      p.prm = tryCatch({coef(summary(prm))[2,3]}, error=function(e){
        print(paste0('error of NA'))
        return(NA)}),
      
      
      
      beta = tryCatch(round(coef(summary(fit_un))[2,1],4), error=function(e){print(paste0('error of tissue '))
        return(NA)}),
      
      #StdEr = round(coef(summary(fit_un))[2,2], digits = 4),
      
      p.value_adj_sex_age = coef(summary(fit_adj))[2,4],
      p.prm_adj_sex_age = coef(summary(prm_adj))[2,3],
      beta_adj_sex_age = coef(summary(fit_adj))[2,1],
      StdEr_adj_sex_age = round(coef(summary(fit_adj))[2,2], digits = 4))
    
    
    # bind rows of temporary data frame to the results data frame
    
    df.out <- rbind(df.out, df.lm)
    
    
  }
}



##################################
# Test the Boxplot 
##################################

# 3 SNPs p-value is sig, try to do plot 
# rs1374468_add
# rs2854915_add
# rs11930034_add


# Use the data used for lm() analysis 
# inputdf

# vln box plot for significant gene TPM with 17SNPs from lm result
library(ggplot2)
library(ggsci)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(rstatix)
library(dplyr)
library(scales)
library(patchwork)

# try to do plot 
# use LINC01426 as example 
# i <- "LINC01426"

sig_snp <- c('rs1374468','rs2854915','rs11930034')

l.plot <- inputdf %>% dplyr::select(Sample_Name,FGFR3,Ancestry,rs1374468,rs1063743_add,rs2854915,rs2854915_add,rs11930034,rs11930034_add) 

# ############################ FOR changing GT order ####################################################################
# 
# 
# for (i in 1:nrow(l.sig)){
# 
# # i <- 11
# 
# sp <- l.sig$snp[i]
# gene <- l.sig$variable[i]
# 
# l.plot_order <- l.plot %>%
#   filter(!is.na(l.plot[[sp]]) & l.plot[[sp]] != "")
# 
# s1 <- c(l.sig$geno_1[i],l.sig$geno_2[i])
# 
# if(sp == "rs7279257"){
#   l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = c('TT','AT','AA'))
#   } else {l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = s1)}
# 
# 
# 
# p <-l.plot_order %>%
#   ggplot(aes(x = .data[[sp]], y = .data[[gene]], fill = .data[[sp]])) +
#   geom_violin(width = 1, lwd = 0.2, colour = "black",
#               show.legend = NA, inherit.aes = TRUE) +
#   # geom_violin(trim = FALSE) +
#   scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),drop=F) + ## modified
#   
#   geom_boxplot(
#     width=0.15,lwd = 0.5,fill="white", position = position_dodge(width = 1),
#     outlier.shape=NA) + 
#   
#   geom_point(shape = 1, size = 0.75, colour = "black",
#              position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) +
#   # Add the mean point
#   stat_summary(fun = mean, geom = "point",
#                position = position_dodge(width = 1),
#                shape = 19, size = 1.5, colour = "red") + theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.y = element_text(color = "black",size = 10),
#         axis.text.x = element_text(color = "black",size = 10), 
#         axis.title.y =element_blank(),legend.position = "none") + 
#   ggtitle(paste(sp,'with expression of',gene,'\n p-value:',round(l.sig$p.value_adj_sex_age_ebv[i],3),', beta:',round(l.sig$beta_adj_sex_age_ebv[i],3))) 
# 
# 
# 
# print(p)
# 
# }
# 
# 
# 

########### plot rs11145 and rs224  ########### 

# cahnge order of rs111457485
sig_snp
plots <- list()
gene <- "FGFR3"
sp <- sig_snp[3]

for (i in gene){
  
  # i <- "FGFR3"
  counts <- l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "") %>%  group_by(.data[[sp]]) %>% summarise(n=n())
  x_labels <- paste(counts[[sp]], "\n(n=", counts$n, ")", sep = "")
  # y_max <- max(l.plot[[i]], na.rm = TRUE) + 0.25
  
  
  p <-l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "") %>% 
    ggplot(aes(x = .data[[sp]], y = .data[[i]], fill = .data[[sp]])) +
    geom_boxplot(
      width=0.5,lwd = 0.5, 
      outlier.shape=NA) + 
    
    scale_fill_manual(values = c("#F38491", "#6ECFF6", "#4197EC"),drop=F) +
    
    geom_point(shape = 1, size = 1.5, colour = "black",
               position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1)) +
    # Add the mean point
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 1),
                 shape = 19, size = 1.5, colour = "red") + theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color = "black",size = 10),
          axis.text.x = element_text(color = "black",size = 8,angle = 45, hjust = 1), 
          axis.title.y =element_blank(),legend.position = "none",
          plot.title = element_text(size = 10, face = "plain", color = "black")) + scale_x_discrete(labels = x_labels) + #ylim(0, 45) +
    ggtitle(paste0("Expression of ",i," vs ",sp))
  
  
  #print(p)
  plots[[i]] <- p
  
}

################################################################################
# CUT OFF
#################################################################################



######################################## for rs224, need to change the GT 
################################################################################
################################################################################

sp <- "rs2242780"
sp <- "rs7279257_md"
gene <- siggene[2:3]
# i <- gene[1]
plots <- list()

l.plot_order <- l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "")
l.plot_order$rs7279257_md <- car::recode(l.plot_order$rs7279257, " 'TT'='AT'   ")
#l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = c("GG","CG"),ordered = T)
l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = c("CC","CT"),ordered = T)

for (i in gene){
  counts <- l.plot_order %>%  group_by(.data[[sp]]) %>% summarise(n=n())
  x_labels <- paste(counts[[sp]], "\n(n=", counts$n, ")", sep = "")
  # y_max <- max(l.plot[[i]], na.rm = TRUE) + 0.25
  
  
  p <-l.plot_order %>% 
    ggplot(aes(x = .data[[sp]], y = .data[[i]], fill = .data[[sp]])) +
    geom_boxplot(
      width=0.5,lwd = 0.5, 
      outlier.shape=NA) + 
    
    scale_fill_manual(values = c("#F38491", "#6ECFF6", "#4197EC"),drop=F) +
    
    geom_point(shape = 1, size = 1.5, colour = "black",
               position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1)) +
    # Add the mean point
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 1),
                 shape = 19, size = 1.5, colour = "red") + theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color = "black",size = 10),
          axis.text.x = element_text(color = "black",size = 8,angle = 45, hjust = 1), 
          axis.title.y =element_blank(),legend.position = "none",
          plot.title = element_text(size = 10, face = "plain", color = "black")) + scale_x_discrete(labels = x_labels) +
    ggtitle(paste0("Expression of ",i," vs ",sp))
  
  
  #print(p)
  #plots[[i]] <- p
  plots[[paste0(i,"__rs5")]] <- p
}

combined_plot <- plots[["LINC01426"]] + plots[["LINC01426_ENST00000420877.1"]] +  plots[["LINC01426__rs5"]] + plots[["LINC01426_ENST00000420877.1__rs5"]] + plot_layout(ncol = 4)
print(combined_plot)

