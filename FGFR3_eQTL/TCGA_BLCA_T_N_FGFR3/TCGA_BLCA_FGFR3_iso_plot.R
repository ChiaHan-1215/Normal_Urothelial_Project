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

BLCA_N_fgfr_sex <- left_join(BLCA_N_fgfr,TCGA_sex,by="PID")
BLCA_N_fgfr_sex <- BLCA_N_fgfr_sex[,c(1,13,2:12)]


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


blca_N_long_s %>%
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
  scale_fill_manual(values = c("MALE" = "#A8D5E2", "FEMALE" = "#F2BAC9")) +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
        legend.position = "top")

