# Goal: Do SummaryStat of GTEx v10 bladder tissue FGFR3 isoform TPM 


### Load VcfR 
library(vcfR)
library(biomaRt)
library(dplyr)
library(tidyr)



# load the sample isofom count
iso <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/GTEx_v10_FGFR3_iso_TPM.txt')
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


##########################################################################################

# get TPM
#tmp.i_tpm <- ls_tpm[[i]]
#head(tmp.i_tpm[,1:10])
#tmp.i_tpm <- tmp.i_tpm[ which(tmp.i_tpm$Name %in% tmp.i_tmm$Name), ]
#tmp.i_tpm <- tmp.i_tpm %>% mutate(across(3:ncol(tmp.i_tpm), ~ as.numeric(as.character(.x))))

mani_bl_stat <- t(mani_bladder_with_iso) %>% as.data.frame()
mani_bl_stat <- mani_bl_stat[-1,]
mani_bl_stat <- mani_bl_stat %>% mutate(GENEID=rownames(mani_bl_stat),.before = 1)
mani_bl_stat <- mani_bl_stat %>% mutate(across(2:ncol(mani_bl_stat), ~ as.numeric(as.character(.x))))

str(mani_bl_stat)

tmp.i_tpm <- mani_bl_stat


df_stats.tpm <- tmp.i_tpm %>%
  rowwise() %>%
  mutate(
    # Specify the range of your actual data columns
    min = min(c_across(2:ncol(tmp.i_tpm)), na.rm = TRUE),
    max = max(c_across(2:ncol(tmp.i_tpm)), na.rm = TRUE),
    mean = mean(c_across(2:ncol(tmp.i_tpm)), na.rm = TRUE),
    median = median(c_across(2:ncol(tmp.i_tpm)), na.rm = TRUE),
    sd = sd(c_across(2:ncol(tmp.i_tpm)), na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  select(min, max, mean, median, sd)

# add gene ID, symbol, and method tag
df_stats.tpm <- bind_cols(tmp.i_tpm[,1], df_stats.tpm)
names(df_stats.tpm)[2:ncol(df_stats.tpm)] <- paste0(names(df_stats.tpm)[2:ncol(df_stats.tpm)], "_TPM")
names(df_stats.tpm)[1] <- "FGFR3_isoform"
