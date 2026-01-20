
### Load VcfR 
library(vcfR)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

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



# load the sample TPM count
iso <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/GTEx_v10_FGFR3_iso_TPM.txt')

# load the sample isofom count
# iso <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/FGFR3_iso.txt')

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
mani_bladder_with_iso <- mani_bladder_with_iso[,c(-12)]

#### Vln Box plot for isoform expression between baldder vs smooth muscle ####

# datasets
mani_tot <- mani[,c(2,12,14)]
names(mani_tot)[1] <- "Sample_ID"
iso_p <- left_join(iso,mani_tot,by='Sample_ID')
iso_p <- iso_p %>% select("Sample_ID","SMTS","SMTSD",, everything()) 

############# Do summaryStat of ALL ###################


# use for loop
tissue <- iso_p$SMTS %>% unique() %>% sort()

Summ <- data.frame()

for (j in tissue){
  # j <- "Bladder"
  # j <- "Ovary"
  
  df1 <- iso_p %>% filter(SMTS == j)
  df1 <- t(df1) %>% as.data.frame()
  names(df1) <- df1[1,]
  df1 <- df1[-c(1:3),]
  df1 <- df1 %>% mutate(GENEID=rownames(df1),.before = 1)
  df1 <- df1 %>% mutate(across(2:ncol(df1), ~ as.numeric(as.character(.x))))
  
  df_stats.tpm <- df1 %>%
    rowwise() %>%
    mutate(
      # Specify the range of your actual data columns
      min = min(c_across(2:ncol(df1)), na.rm = TRUE),
      max = max(c_across(2:ncol(df1)), na.rm = TRUE),
      mean = mean(c_across(2:ncol(df1)), na.rm = TRUE),
      median = median(c_across(2:ncol(df1)), na.rm = TRUE),
      sd = sd(c_across(2:ncol(df1)), na.rm = TRUE)
    ) %>%
    ungroup() %>% 
    select(min, max, mean, median, sd)
  
  # add gene ID, symbol, and method tag
  df_stats.tpm <- bind_cols(df1[,1], df_stats.tpm)
  names(df_stats.tpm)[2:ncol(df_stats.tpm)] <- paste0(names(df_stats.tpm)[2:ncol(df_stats.tpm)], "_TPM")
  names(df_stats.tpm)[1] <- "FGFR3_isoform"
  df_stats.tpm <- df_stats.tpm %>% mutate(Tissue=j,.before = 1)
  
  Summ <- rbind(df_stats.tpm,Summ)
  rm(df_stats.tpm)
  rm(df1)
  
  
}

############# ############# ############# ############# 
############# ############# ############# ############# 
############# ############# ############# ############# 

# Now plot for bladder and muscle
iso_p_bm <- iso_p %>% filter(SMTS=="Bladder"|SMTS=="Muscle")

#
iso_long <- iso_p_bm %>%
  pivot_longer(
    cols = -c(Sample_ID, SMTS, SMTSD),
    names_to = "isoform",
    values_to = "tpm"
  )

# counts per tissue for axis labels
counts <- iso_long %>%
  distinct(Sample_ID, SMTS) %>%
  count(SMTS, name = "n")

iso_long2 <- iso_long %>%
  left_join(counts, by = "SMTS") %>%
  mutate(myaxis = paste0(SMTS, ", n=", n))
# 
# 
# ggplot(iso_long2, aes(x = factor(myaxis, levels = rev(unique(myaxis))), y = tpm)) +
#   geom_violin(width = 0.8, linewidth = 0.15, colour = "black", 
#               fill = "#F2C5B3", trim = F, alpha = 0.6) +
#   geom_jitter(width = 0.15, size = 0.1, alpha = 0.3, color = "grey30") +
#   stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000") +
#   facet_wrap(~ isoform, nrow = 5, ncol = 2, scales = "free_y") +
#   theme_classic() +
#   theme(axis.title = element_blank(),
#         axis.text.y = element_text(color = "black", size = 10),
#         legend.position = "none")
# 
# 
# 


# Muscle plot
p_muscle <- iso_long2 %>%
  filter(SMTS == "Muscle") %>%
  ggplot(aes(x = isoform, y = log2(tpm+1))) +  # Changed x to isoform
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
  ggtitle("Muscle (n=818)")


p_bladder <- iso_long2 %>%
  filter(SMTS == "Bladder") %>%
  ggplot(aes(x = isoform, y = log2(tpm+1))) +  # Changed x to isoform
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
  ggtitle("Bladder (n=77)")  # Add sample size manually


p_bladder  | p_muscle
