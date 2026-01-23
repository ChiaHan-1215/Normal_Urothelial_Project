
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

# Now plot for bladder 
iso_p_bm <- iso_p %>% filter(SMTS=="Bladder")

# Load the GTEx gender info
gtex.gender <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/Metadata_Files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt')
names(gtex.gender)[1] <- "GTEx_ID"
gtex.gender <- gtex.gender[,c(1,2)]

iso_p_bm <- iso_p_bm %>% mutate(GTEx_ID=sub("^([^-]+-[^-]+).*$", "\\1", iso_p_bm$Sample_ID),.before = 1)

iso_p_bm <- inner_join(gtex.gender,iso_p_bm,by="GTEx_ID")



#
iso_long <- iso_p_bm %>%
  pivot_longer(
    cols = -c(GTEx_ID,SEX,Sample_ID, SMTS, SMTSD),
    names_to = "isoform",
    values_to = "tpm"
  )

# counts per tissue for axis labels
counts <- iso_long %>%
  distinct(Sample_ID, SEX) %>%
  count(SEX, name = "n")

iso_long2 <- iso_long %>%
  left_join(counts, by = "SEX") %>%
  mutate(myaxis = paste0(SEX, ", n=", n))

# plotting 

iso_long2 %>%
  ggplot(aes(x = isoform, y = log2(tpm + 1), fill = factor(SEX))) +  # 1. Added fill here
  geom_boxplot(width = 0.6, 
               linewidth = 0.5, 
               colour = "black", 
               alpha = 0.6, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) + # 2. Optional: adjust spacing
  geom_jitter(aes(group = factor(SEX)),                  # 3. Group dots by sex
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
  scale_fill_manual(values = c("1" = "#A8D5E2", "2" = "#F2BAC9"), 
                    labels = c("Male", "Female"),
                    name = "Sex") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "top")

