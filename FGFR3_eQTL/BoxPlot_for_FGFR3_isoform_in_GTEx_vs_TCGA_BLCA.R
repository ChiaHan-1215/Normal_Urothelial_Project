### Gaol: plotting FGFR3 isoform expression for TCGA BLCA T vs N

library(vcfR)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/Original source data Expression/UCSC_Toil_RNAseq_Recompute/transcript expression RNAseq_hg38/byCancerType/')

# load the sample TPM count
BLCA_N <- read.delim('TCGA_BLCA_IsoformExpression.Profile_OnlyNormal.txt')
BLCA_T <- read.delim('TCGA_BLCA_IsoformExpression.Profile_OnlyTumor.txt')

# subset FGFR3
BLCA_N_fgfr <- BLCA_N %>% filter(gene == 'FGFR3')
BLCA_T_fgfr <- BLCA_T %>% filter(gene == 'FGFR3')

# As the TPM value is transfomed, need to reconvert to TPM
#log2(tpm+0.001)

BLCA_N_fgfr[, 7:ncol(BLCA_N_fgfr)] <- apply(BLCA_N_fgfr[, 7:ncol(BLCA_N_fgfr)], 
                                            MARGIN = c(1, 2), 
                                            FUN = function(x) round(2^x - 0.001, 2))

BLCA_T_fgfr[, 7:ncol(BLCA_T_fgfr)] <- apply(BLCA_T_fgfr[, 7:ncol(BLCA_T_fgfr)], 
                                            MARGIN = c(1, 2), 
                                            FUN = function(x) round(2^x - 0.001, 2))


# For Normal
BLCA_N_fgfr <- data.frame(t(BLCA_N_fgfr))
names(BLCA_N_fgfr) <-BLCA_N_fgfr[1,] 
BLCA_N_fgfr <- BLCA_N_fgfr[-c(1:6),]

BLCA_N_fgfr <- BLCA_N_fgfr %>%
  mutate(across(everything(), as.numeric))

BLCA_N_fgfr <- BLCA_N_fgfr %>% mutate(PID=rownames(BLCA_N_fgfr),.before = 1)
BLCA_N_fgfr$total_FGFR3 <- rowSums(BLCA_N_fgfr[,c(2:11)])
BLCA_N_fgfr$PID <- gsub('\\.','-',BLCA_N_fgfr$PID)


# For tumor 
BLCA_T_fgfr <- data.frame(t(BLCA_T_fgfr))
names(BLCA_T_fgfr) <-BLCA_T_fgfr[1,] 
BLCA_T_fgfr <- BLCA_T_fgfr[-c(1:6),]

BLCA_T_fgfr <- BLCA_T_fgfr %>%
  mutate(across(everything(), as.numeric))

BLCA_T_fgfr <- BLCA_T_fgfr %>% mutate(PID=rownames(BLCA_T_fgfr),.before = 1)
BLCA_T_fgfr$total_FGFR3 <- rowSums(BLCA_T_fgfr[,c(2:11)])
BLCA_T_fgfr$PID <- gsub('\\.','-',BLCA_T_fgfr$PID)


#### Vln Box plot for isoform expression between baldder vs smooth muscle ####

blca_T_long <- BLCA_T_fgfr %>%
  pivot_longer(
    cols = -c(PID,total_FGFR3),
    names_to = "isoform",
    values_to = "tpm"
  )


blca_N_long <- BLCA_N_fgfr %>%
  pivot_longer(
    cols = -c(PID,total_FGFR3),
    names_to = "isoform",
    values_to = "tpm"
  )


# Muscle plot
p_blca_Tumor <- blca_T_long %>%
  ggplot(aes(x = isoform, y = log2(tpm+1))) +  # Changed x to isoform
  geom_boxplot(width = 0.4, linewidth = 0.5, colour = "black", 
               fill = "#A8D5E2", alpha = 0.6,outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.1, alpha = 0.5, color = "grey30") +
  stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000") +

  # Removed facet_wrap
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "none") + ylim(c(0,12)) +
  ggtitle("Tumor (n=407)")


p_blca_Normal <- blca_N_long %>%
  ggplot(aes(x = isoform, y = log2(tpm+1))) +  # Changed x to isoform
  geom_boxplot(width = 0.4, linewidth = 0.5, colour = "black", 
               fill = "#A8D5E2", alpha = 0.6,outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 0.1, alpha = 0.5, color = "grey30") +
  stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000") +
  
  # Removed facet_wrap
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "none") + ylim(c(0,12)) +
  ggtitle("Normal (n=19)")


p_blca_Normal  | p_blca_Tumor




####Adding Sample SEX info


TCGA_sex <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/TCGA_masterFiles/TCGA master table 2017.xlsx',sheet = 2)

TCGA_sex <- TCGA_sex[,c(1,11)]
names(TCGA_sex)[1] <- "PID"

# Now we have TN data
BLCA_T_fgfr
BLCA_N_fgfr


BLCA_T_fgfr_sex <- left_join(BLCA_T_fgfr,TCGA_sex,by="PID")
BLCA_T_fgfr_sex <- BLCA_T_fgfr_sex[,c(1,13,2:12)]
BLCA_T_fgfr_sex$gender <- factor(BLCA_T_fgfr_sex$gender,levels = c("MALE","FEMALE"))


BLCA_N_fgfr_sex <- left_join(BLCA_N_fgfr,TCGA_sex,by="PID")
BLCA_N_fgfr_sex <- BLCA_N_fgfr_sex[,c(1,13,2:12)]
BLCA_N_fgfr_sex$gender <- factor(BLCA_N_fgfr_sex$gender,levels = c("MALE","FEMALE"))




# PLOT 

blca_T_long_s <- BLCA_T_fgfr_sex %>%
  pivot_longer(
    cols = -c(PID,total_FGFR3,gender),
    names_to = "isoform",
    values_to = "tpm"
  )


blca_N_long_s <- BLCA_N_fgfr_sex %>%
  pivot_longer(
    cols = -c(PID,total_FGFR3,gender),
    names_to = "isoform",
    values_to = "tpm"
  )





####**************************************************************************************

# choose isoforms of interest
keep_isoforms <- c(
  "ENST00000340107.8",
  "ENST00000481110.6",
  "ENST00000352904.5"
)

blca_N_long_s_filt <- blca_N_long_s %>%
  filter(isoform %in% keep_isoforms)  %>%
  mutate(isoform = factor(isoform, levels = keep_isoforms))

blca_T_long_s_filt <- blca_T_long_s %>%
  filter(isoform %in% keep_isoforms) %>%
  mutate(isoform = factor(isoform, levels = keep_isoforms))

####**************************************************************************************


t_normal <- blca_N_long_s_filt %>%
  ggplot(aes(x = isoform, y = log2(tpm + 1), fill = gender)) +  # 1. Added fill here
  geom_boxplot(width = 0.6, 
               linewidth = 0.5, 
               colour = "black", 
               alpha = 0.6, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) + # 2. Optional: adjust spacing
  geom_jitter(aes(group = gender),                  # 3. Group dots by sex
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
  
  stat_compare_means(method = "t.test", 
                     label = "p.format") +
  
  
  
  
  scale_fill_manual(values = c("MALE" = "#A8D5E2", "FEMALE" = "#F2BAC9"),labels = c("M,n=10", "F,n=9")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "right") + ylim(c(-0.2,10.5))+ ggtitle("BLCA normal")




t_Tumor <- blca_T_long_s_filt %>%
  ggplot(aes(x = isoform, y = log2(tpm + 1), fill = gender)) +  # 1. Added fill here
  geom_boxplot(width = 0.6, 
               linewidth = 0.5, 
               colour = "black", 
               alpha = 0.6, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) + # 2. Optional: adjust spacing
  geom_jitter(aes(group = gender),                  # 3. Group dots by sex
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
  
  stat_compare_means(method = "t.test", 
                     label = "p.format") +
  
  
  scale_fill_manual(values = c("MALE" = "#A8D5E2", "FEMALE" = "#F2BAC9"),labels = c("M,n=301", "F,n=106")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "right") + ylim(c(-0.2,10.5))+ ggtitle("BLCA tumor")




#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################

# FOR GTEX 

#####################################################################################################################################
#####################################################################################################################################
#####################################################################################################################################


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

#******

GTEx_keep_isoforms <- c('ENST00000340107.8','ENST00000481110.7','ENST00000352904.6')

iso_long2_filt <- iso_long2 %>%
  filter(isoform %in% GTEx_keep_isoforms) %>%
  mutate(isoform = factor(isoform, levels = GTEx_keep_isoforms))

#******

library(ggpubr)

g1 <- iso_long2_filt %>%
  ggplot(aes(x = isoform, y = log2(tpm + 1), fill = factor(SEX))) +  # 1. Added fill here
  geom_boxplot(width = 0.6, 
               linewidth = 0.5, 
               colour = "black", 
               alpha = 0.6, 
               outlier.shape = NA,
               position = position_dodge(width = 0.8)) + 
  geom_jitter(aes(group = factor(SEX)),                  
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), 
              size = 0.1, 
              alpha = 0.3, 
              color = "grey30") +
  
  stat_summary(
               fun = median, 
               geom = "point", 
               shape = 20, 
               size = 2, 
               color = "#FF0000",
               position = position_dodge(width = 0.8)) +
  
  stat_compare_means(method = "t.test", 
                     label = "p.format") +

  
  scale_fill_manual(values = c("1" = "#A8D5E2", "2" = "#F2BAC9"), 
                    labels = c("M,n=48", "F,n=29"),
                    name = "gender") +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "right")+ ylim(c(-0.2,10.5)) + ggtitle("GTEx") 
  





# Save plot 9 * 4.2

t_normal|t_Tumor
g1|t_normal
g1|t_Tumor

## check from dataset 
# Creates a table of p-values for every isoform
check_stats <- iso_long2_filt %>%
  group_by(isoform) %>%
  summarise(
    t_test_p = t.test(log2(tpm + 1) ~ SEX)$p.value,
    wilcox_p = wilcox.test(tpm ~ SEX)$p.value
  )

print(check_stats)


