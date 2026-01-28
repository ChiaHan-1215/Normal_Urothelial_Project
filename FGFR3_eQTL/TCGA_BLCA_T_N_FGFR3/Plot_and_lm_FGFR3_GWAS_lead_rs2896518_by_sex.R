### Gaol: plotting FGFR3 isoform expression for TCGA BLCA T vs N

library(vcfR)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(lmPerm)


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




GT_data1 <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/Data_for_mapBC_project/mapBC_genotypes_TCGA_BLCA_72markers.csv')
GT_data2 <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/Data_for_mapBC_project/mapBC_genotypes_rs2236786_TCGA_BLCA_Tumor_Blood.csv')

# GWAS lead rs2896518
GT_data1 <- GT_data1[,c("ID_person","rs2896518")]
names(GT_data1)[1] <- "PID"



# merge data
BLCA_N_fgfr <- left_join(BLCA_N_fgfr,GT_data,by="PID")
BLCA_T_fgfr <- left_join(BLCA_T_fgfr,GT_data,by="PID")


# add 0,1,2

# Now add 0,1,2 to the GT 

add_01_GT <- function(inputdata){
  
  missing_gt <- c(".", "", "NA")
  
  for (i in grep('rs',names(inputdata),value = T)) {
    #i <- "rs372189430"
    n1 <- paste0(i, "_add")
    
    # --- 0) Normalize missing genotype codes to NA ---
    x <- inputdata[[i]]
    x[x %in% missing_gt] <- NA
    inputdata[[i]] <- x  # write back so downstream uses cleaned data
    
    
    # 1. Create frequency table and ensure it's a clean data frame
    tp1 <- table(inputdata[[i]], useNA = "no")
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
      inputdata[[n1]] <- car::recode(inputdata[[i]], recode_str)
      
    } else if (nrow(SAMEGT) == 1) {
      # Case with only 1 homozygote (recode to 0) and potentially 1 heterozygote (recode to 1)
      if (nrow(tp_left) >= 1) {
        recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0; '", tp_left$Genotype[1], "'=1")
      } else {
        recode_str <- paste0("'", SAMEGT$Genotype[1], "'=0")
      }
      inputdata[[n1]] <- car::recode(inputdata[[i]], recode_str)
    }
  }
  return(inputdata)
}


BLCA_N_fgfr_GT <- add_01_GT(inputdata =BLCA_N_fgfr )
BLCA_T_fgfr_GT <- add_01_GT(inputdata =BLCA_T_fgfr )


BLCA_N_fgfr_GT <- BLCA_N_fgfr_GT[,c(1,13,14,2:11)]
BLCA_T_fgfr_GT <- BLCA_T_fgfr_GT[,c(1,13,14,2:11)]



# 
# #### Vln Box plot ####
# 
# blca_T_long <- BLCA_T_fgfr %>%
#   pivot_longer(
#     cols = -c(PID,total_FGFR3),
#     names_to = "isoform",
#     values_to = "tpm"
#   )
# 
# 
# blca_N_long <- BLCA_N_fgfr %>%
#   pivot_longer(
#     cols = -c(PID,total_FGFR3),
#     names_to = "isoform",
#     values_to = "tpm"
#   )
# 
# 
# # Muscle plot
# p_blca_Tumor <- blca_T_long %>%
#   ggplot(aes(x = isoform, y = log2(tpm+1))) +  # Changed x to isoform
#   geom_boxplot(width = 0.4, linewidth = 0.5, colour = "black", 
#                fill = "#A8D5E2", alpha = 0.6,outlier.shape = NA) +
#   geom_jitter(width = 0.15, size = 0.1, alpha = 0.5, color = "grey30") +
#   stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000") +
#   # Removed facet_wrap
#   theme_classic() +
#   theme(axis.title = element_blank(),
#         axis.text.y = element_text(color = "black", size = 10),
#         axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
#         legend.position = "none") + ylim(c(0,12)) +
#   ggtitle("Tumor (n=407)")
# 
# 
# p_blca_Normal <- blca_N_long %>%
#   ggplot(aes(x = isoform, y = log2(tpm+1))) +  # Changed x to isoform
#   geom_boxplot(width = 0.4, linewidth = 0.5, colour = "black",
#                fill = "#A8D5E2", alpha = 0.6,outlier.shape = NA) +
#   geom_jitter(width = 0.15, size = 0.1, alpha = 0.5, color = "grey30") +
#   stat_summary(fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000") +
#   # Removed facet_wrap
#   theme_classic() +
#   theme(axis.title = element_blank(),
#         axis.text.y = element_text(color = "black", size = 10),
#         axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
#         legend.position = "none") + ylim(c(0,12)) +
#   ggtitle("Normal (n=19)")
# 
# 
# p_blca_Normal  | p_blca_Tumor
# 
# 


####Adding Sample SEX info


TCGA_sex <- readxl::read_xlsx('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/TCGA_data/TCGA_masterFiles/TCGA master table 2017.xlsx',sheet = 2)

TCGA_sex <- TCGA_sex[,c(1,11)]
names(TCGA_sex)[1] <- "PID"

# Now we have TN data
BLCA_T_fgfr_GT
BLCA_N_fgfr_GT



BLCA_T_fgfr_sex <- left_join(BLCA_T_fgfr_GT,TCGA_sex,by="PID")
BLCA_T_fgfr_sex <- BLCA_T_fgfr_sex[,c(1,14,2:13)]

BLCA_N_fgfr_sex <- left_join(BLCA_N_fgfr_GT,TCGA_sex,by="PID")
BLCA_N_fgfr_sex <- BLCA_N_fgfr_sex[,c(1,14,2:13)]

BLCA_T_fgfr_sex$gender <- factor(BLCA_T_fgfr_sex$gender,levels = c("MALE","FEMALE"))
BLCA_N_fgfr_sex$gender <- factor(BLCA_N_fgfr_sex$gender,levels = c("MALE","FEMALE"))

# choose isoforms of interest
keep_isoforms <- c(
  "ENST00000340107.8",
  "ENST00000481110.6",
  "ENST00000352904.5"
)

blca_N_filt <- BLCA_N_fgfr_sex %>% select(names(BLCA_N_fgfr_sex)[1:4],all_of(keep_isoforms))
blca_T_filt <- BLCA_T_fgfr_sex %>% select(names(BLCA_T_fgfr_sex)[1:4],all_of(keep_isoforms))



# PLOT 

blca_N_long <- blca_N_filt %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "isoform",
    values_to = "tpm"
  ) %>%
  mutate(
    gender = factor(gender, levels = c("MALE","FEMALE")),
    rs2896518 = factor(rs2896518)  # A/G, G/G, etc.
  )


blca_T_long <- blca_T_filt %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "isoform",
    values_to = "tpm"
  ) %>%
  mutate(
    gender = factor(gender, levels = c("MALE","FEMALE")),
    rs2896518 = factor(rs2896518)  # A/G, G/G, etc.
  )



p_blca_sex_gt_N <- blca_N_long %>%
  ggplot(aes(x = rs2896518, y = log2(tpm + 1), fill = gender)) +
  geom_boxplot(
    width = 0.4, linewidth = 0.5, colour = "black",
    alpha = 0.6, outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(group = gender),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 0.1, alpha = 0.5, color = "grey30"
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000",
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = c("MALE" = "#A8D5E2", "FEMALE" = "#F2BAC9")) +
  facet_wrap(~ isoform, scales = "free_y") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  ggtitle("BLCA Normal: rs2896518 genotype vs FGFR3 isoform expression")


p_blca_sex_gt_T <- blca_T_long %>%
  ggplot(aes(x = rs2896518, y = log2(tpm + 1), fill = gender)) +
  geom_boxplot(
    width = 0.4, linewidth = 0.5, colour = "black",
    alpha = 0.6, outlier.shape = NA,
    position = position_dodge(width = 0.8)
  ) +
  geom_jitter(
    aes(group = gender),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8),
    size = 0.1, alpha = 0.5, color = "grey30"
  ) +
  stat_summary(
    fun = median, geom = "point", shape = 20, size = 2, color = "#FF0000",
    position = position_dodge(width = 0.8)
  ) +
  scale_fill_manual(values = c("MALE" = "#A8D5E2", "FEMALE" = "#F2BAC9")) +
  facet_wrap(~ isoform, scales = "free_y") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 9, angle = 45, hjust = 1),
    legend.position = "right"
  ) +
  ggtitle("BLCA Tumor: rs2896518 genotype vs FGFR3 isoform expression")




p_blca_sex_gt_N|p_blca_sex_gt_T



### Do lm() 

inputdf <- blca_N_filt
inputdf <- blca_T_filt

df.out <- data.frame()
df_count.summary <- data.frame()

for (i in grep("_add",names(inputdf),value = T)){

  # i <- "rs2896518_add"

  
  
  function.names <- c("max", "min", "mean", "median", "sd")
  
  for ( k in function.names ) {
    # k <- function.names[3]
    # the all the TPM in all tissue
    # consider dealing with NA
    
    # & inputdf[[i]] != 0
    # names(inputdf[6:135]
    
    # remove SNPs that are NA and empty
    df.tmp <- inputdf[,c(names(inputdf[5:7]), gsub("_add", "", i))] %>%
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
  
  for(j in names(inputdf[5:7])){
    
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
    fmla_adj <- as.formula(paste0(j, "~" , i, " + gender "))
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





##### CUT OF 


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
                                 

# 
# # cahnge order of rs111457485
# sig_snp
# plots <- list()
# gene <- "FGFR3"
# sp <- sig_snp[3]
# 
# for (i in gene){
#   
#   # i <- "FGFR3"
#   counts <- l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "") %>%  group_by(.data[[sp]]) %>% summarise(n=n())
#   x_labels <- paste(counts[[sp]], "\n(n=", counts$n, ")", sep = "")
#   # y_max <- max(l.plot[[i]], na.rm = TRUE) + 0.25
#   
#   
#   p <-l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "") %>% 
#     ggplot(aes(x = .data[[sp]], y = .data[[i]], fill = .data[[sp]])) +
#     geom_boxplot(
#       width=0.5,lwd = 0.5, 
#       outlier.shape=NA) + 
#     
#     scale_fill_manual(values = c("#F38491", "#6ECFF6", "#4197EC"),drop=F) +
#     
#     geom_point(shape = 1, size = 1.5, colour = "black",
#                position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1)) +
#     # Add the mean point
#     stat_summary(fun = mean, geom = "point",
#                  position = position_dodge(width = 1),
#                  shape = 19, size = 1.5, colour = "red") + theme_classic() + 
#     theme(axis.title.x=element_blank(),
#           axis.text.y = element_text(color = "black",size = 10),
#           axis.text.x = element_text(color = "black",size = 8,angle = 45, hjust = 1), 
#           axis.title.y =element_blank(),legend.position = "none",
#           plot.title = element_text(size = 10, face = "plain", color = "black")) + scale_x_discrete(labels = x_labels) + #ylim(0, 45) +
#     ggtitle(paste0("Expression of ",i," vs ",sp))
#   
#   
#   #print(p)
#   plots[[i]] <- p
