### Goal: Analysis FGFR3 isofroma eQTL using GTEx v10 and Upcoming v11 in bladder tissue data

### Data Need: 
  - As the newer version has around 77 bladder tissue, check those new sample to see if it have GT info
  - The GTEx v10 data are located in T-drive: `/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10`
      - Insdie the folder, the `/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/QTL/eQTL_covariates/Bladder.v10.covariates.txt` maybe useful for adjusting lm() analysis 
   
  - WGS vcf file are in GTEx folder, check to see if the ID match


### Script:

```R

### Load VcfR 
library(vcfR)
library(biomaRt)

# ==============================================================================
# 1) LOAD VCF AND EXTRACT DATA
# ==============================================================================

# load cov

cov <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/QTL/eQTL_covariates/Bladder.v10.covariates.txt')
cov <- t(cov) %>% as.data.frame()
names(cov) <- cov[1,]
cov <- cov[-1,]


cov$Sample_ID <- rownames(cov) %>% gsub('\\.','-',.)


vcf.file.tmp <- read.vcfR('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/GTEx_WGS_subset_FGFR3_TACC3.vcf.gz')

vcf_genotypes.tmp <- extract.gt(
  vcf.file.tmp,
  element = "GT",
  mask = FALSE,
  as.numeric = FALSE,
  return.alleles = TRUE,
  IDtoRowNames = TRUE,
  extract = TRUE,
  convertNA = TRUE
)

geno_df.tmp <- as.data.frame(t(vcf_genotypes.tmp))
# use match function, but %in% also work well.
# %in% give T/F logical output while match gives vector position 
# names(geno_df.tmp) <- snp.list$rsid[ match(gsub("^chr5_|_[A-Z]_.+", "", names(geno_df.tmp)), gsub("chr5:|:[A-Z]", "", snp.list$SNP.pick_hg38)) ]
geno_df.tmp$ID_person <- row.names(geno_df.tmp)
geno_df.tmp$GTEx_ID <- sub("-[^-]+$", "", geno_df.tmp$ID_person)
geno_df.tmp <- geno_df.tmp[,c("ID_person","GTEx_ID",setdiff(names(geno_df.tmp),c("ID_person","GTEx_ID")))]

# remove GTEx-xxx-xxxx-xx to GTEx-xxxx
geno_df.tmp$GTEx_ID <- sub("^([^-]+-[^-]+).*$", "\\1", geno_df.tmp$ID_person)
geno_df.tmp <- geno_df.tmp[,-1]
# md the SNP column names 
# strip the trailing _123
nm <- gsub("_[0-9]+$", "", names(geno_df.tmp))
# make any duplicates unique: x, x_1, x_2, ...
names(geno_df.tmp) <- make.unique(nm, sep = "_")
geno_df.tmp <- data.frame(lapply(geno_df.tmp, function(x){ gsub("\\/", " ", x)}))



# load the sample isofom count
iso <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/FGFR3_iso.txt')
iso <- t(iso) %>% as.data.frame()
names(iso) <- iso[1,]
iso <- iso[-c(1,2),]
iso[] <- lapply(iso, as.numeric)
iso$Sample_ID <- rownames(iso) %>% gsub('\\.','-',.)
iso <- iso[,c(11,1:10)]

# Get the bladder sample ID 

mani <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/Phenotype_dbGaP/phs000424.v10.pht002743.v10.p2.c1.GTEx_Sample_Attributes.GRU.txt',skip = 10)

mani_bladder <- mani %>% dplyr::filter(SMTS == 'Bladder')
mani_bladder <- mani_bladder[,c(1,2)]
names(mani_bladder)[2] <- "Sample_ID"
mani_bladder_with_iso <- inner_join(iso,mani_bladder,by = 'Sample_ID')
rownames(mani_bladder_with_iso) <- mani_bladder_with_iso$Sample_ID
mani_bladder_with_iso <- mani_bladder_with_iso[,c(-1,-12)]


# Now get normalized count 

library(edgeR)

dge <- DGEList(counts = mani_bladder_with_iso)

# calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# get normalized counts
normalized_counts <- cpm(dge)
head(normalized_counts[1:10,1:10])
dim(normalized_counts)

#

inv_norm_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# apply inverse normal transform to each gene
final_normalized <- t(apply(normalized_counts, 1, inv_norm_transform))
final_normalized <- as.data.frame(final_normalized)
head(final_normalized[1:5,1:10])



final_normalized <-final_normalized %>% mutate(Sample_ID=rownames(final_normalized),.before = 1)
final_normalized <-final_normalized %>% mutate(GTEx_ID=sub("^(([^-]+-[^-]+)).*$", "\\1", final_normalized$Sample_ID),.before = 2)


# merge data 
Mg_data <- left_join(final_normalized,geno_df.tmp,by='GTEx_ID')

# Use portion as test, but can be whole dataset after 
Mg_data_s <- Mg_data[,c(1:8,grep('rs',names(Mg_data)))]

# Now add 0,1,2 to the GT 
for (i in grep('rs',names(Mg_data),value = T)) {
  #i <- "rs3135861"
  n1 <- paste0(i, "_add")
  
  # 1. Create frequency table and ensure it's a clean data frame
  tp1 <- table(Mg_data[[i]], useNA = "no")
  if (length(tp1) == 0) next # Skip if column is empty
  
  tp1 <- as.data.frame(tp1, stringsAsFactors = FALSE)
  colnames(tp1) <- c("Genotype", "Freq")
  
  # 2. Split Genotype into Alleles (e.g., "G/G" -> "G", "G")
  tp1 <- tp1 %>%
    separate(Genotype, into = c("A1", "A2"), sep = " ", 
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
    Mg_data_s[[n1]] <- car::recode(Mg_data_s[[i]], recode_str)
    
  } else if (nrow(SAMEGT) == 1) {
    # Case with only 1 homozygote (recode to 0) and potentially 1 heterozygote (recode to 1)
    if (nrow(tp_left) >= 1) {
      recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0; '", tp_left$Genotype[1], "'=1")
    } else {
      recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0")
    }
    Mg_data_s[[n1]] <- car::recode(Mg_data_s[[i]], recode_str)
  }
}

Mg_data_s <- Mg_data_s[,c(names(Mg_data_s)[1:8],setdiff(names(Mg_data_s),names(Mg_data_s)[1:8]) %>% sort())]



## Do lm()






```
