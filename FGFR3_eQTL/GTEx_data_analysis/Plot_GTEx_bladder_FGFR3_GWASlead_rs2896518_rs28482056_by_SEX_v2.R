# Goal : lm anaylsis of GTEx bladder data with GT 
### Load VcfR 
library(vcfR)
library(biomaRt)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lmPerm)
library(patchwork)

# ==============================================================================
# 1) LOAD VCF AND EXTRACT DATA
# ==============================================================================

# load cov, for PC studd

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
mani_bladder_with_iso <- mani_bladder_with_iso[,c(-1,-12)]


# Now get normalized count 
final_normalized <- mani_bladder_with_iso


final_normalized <-final_normalized %>% mutate(Sample_ID=rownames(final_normalized),.before = 1)
final_normalized <-final_normalized %>% mutate(GTEx_ID=sub("^(([^-]+-[^-]+)).*$", "\\1", final_normalized$Sample_ID),.before = 2)


# merge data 
Mg_data <- left_join(final_normalized,geno_df.tmp,by='GTEx_ID')
# merge sample sex age info

GTEx_info <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/GTEx_analysis_v10/Metadata_Files/GTEx_Analysis_2022-06-06_v10_Annotations_GTEx_Analysis_2022-06-06_v10_Annotations_SubjectPhenotypesDS.txt')
names(GTEx_info)[1] <- "GTEx_ID"
GTEx_info <- GTEx_info[,c(1,3,4,5)]

Mg_data <- Mg_data %>%
  left_join(GTEx_info, by = "GTEx_ID") %>%
  dplyr::select(GTEx_ID,
                dplyr::all_of(names(GTEx_info)[names(GTEx_info) != "GTEx_ID"]),
                dplyr::everything())

# write.table('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/GTEx_v10_FGFR3_iso_normalized_with_GT.csv',col.names = T,row.names = F,quote = F,sep = ',')



# SELECT GWAS MERKER FIRST TO TEST

# Use portion as test, but can be whole dataset after 
Mg_data_s <- Mg_data[,c(1:15,grep('rs2896518|rs28482056',names(Mg_data)))]

# Now add 0,1,2 to the GT 

missing_gt <- c(".", "", "NA")

for (i in grep('rs',names(Mg_data_s),value = T)) {
  #i <- "rs372189430"
  n1 <- paste0(i, "_add")
  
  # --- 0) Normalize missing genotype codes to NA ---
  x <- Mg_data_s[[i]]
  x[x %in% missing_gt] <- NA
  Mg_data_s[[i]] <- x  # write back so downstream uses cleaned data
  
  
  # 1. Create frequency table and ensure it's a clean data frame
  tp1 <- table(Mg_data_s[[i]], useNA = "no")
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


GTEx_keep_isoforms <- c('ENST00000340107.8','ENST00000481110.7','ENST00000352904.6')

Mg_data_s <- Mg_data_s %>% dplyr::select(GTEx_ID,SEX,AGE,RACE,rs2896518,rs2896518_add,rs28482056,rs28482056_add,all_of(GTEx_keep_isoforms))


### PLOT

GTEx_GT_iso <- Mg_data_s %>%
  rename(
    PID = GTEx_ID,
    gender = SEX
  ) %>%
  mutate(
    gender = factor(gender, levels = c(1, 2)),
    rs2896518 = factor(rs2896518)
  ) %>%
  pivot_longer(
    cols = starts_with("ENST"),
    names_to = "isoform",
    values_to = "tpm"
  )


#rs28482056 
#rs2896518
snp_plot <- "rs2896518"
#snp_plot <- "rs28482056"

###########################################################################################


Male_plot_data <- GTEx_GT_iso %>%
  # Filter data as before
  filter(gender == 1,
         !is.na(.data[[snp_plot]]),
         .data[[snp_plot]] != "") %>%
  
  # Group by the variables that define each boxplot
  group_by(isoform, .data[[snp_plot]]) %>%
  
  # Create the label using paste() with a newline separator
  # This is the key change that correctly creates a two-line string.
  mutate(genotype_label = paste(.data[[snp_plot]], paste0("n=", n(), ""), sep = "\n")) %>%
  # Ungroup to ensure clean plotting
  ungroup()


p_GTEx_M <- ggplot(Male_plot_data, aes(x = genotype_label,
                                  y = log2(tpm + 1),
                                  fill = gender)) +
  geom_boxplot(
    width = 0.4, linewidth = 0.5,
    colour = "black", alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.8, alpha = 0.5,
    color = "grey30"
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 20, size = 2, color = "#FF0000"
  ) +
  scale_fill_manual(values = c("1" = "#A8D5E2")) +
  facet_wrap(~ isoform, scales = "free_y") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    # For two-line labels, horizontal text is usually best.
    axis.text.x = element_text(color = "black", size = 9, angle = 0, hjust = 1),
    legend.position = "none"
  ) +
  scale_x_discrete(labels = function(x) gsub("\\\\n", "\n", x)) + 
  ggtitle(paste0("Male | ", snp_plot, " genotype vs FGFR3 isoform expression"))

# Display the final plot
print(p_GTEx_M)

#---
  

Female_plot_data <- GTEx_GT_iso %>%
  # Filter data as before
  filter(gender == 2,
         !is.na(.data[[snp_plot]]),
         .data[[snp_plot]] != "") %>%
  
  # Group by the variables that define each boxplot
  group_by(isoform, .data[[snp_plot]]) %>%
  
  # Create the label using paste() with a newline separator
  # This is the key change that correctly creates a two-line string.
  mutate(genotype_label = paste(.data[[snp_plot]], paste0("n=", n(), ""), sep = "\n")) %>%
  # Ungroup to ensure clean plotting
  ungroup()


 
p_GTEx_F <- ggplot(Female_plot_data, aes(x = genotype_label,
                                       y = log2(tpm + 1),
                                       fill = gender)) +
  
  geom_boxplot(
    width = 0.4, linewidth = 0.5,
    colour = "black", alpha = 0.6,
    outlier.shape = NA
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.8, alpha = 0.8,
    color = "grey30"
  ) +
  stat_summary(
    fun = median, geom = "point",
    shape = 20, size = 2, color = "#FF0000"
  ) +
  scale_fill_manual(values = c("2" = "#F2BAC9")) +
  facet_wrap(~ isoform, scales = "free_y") +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(color = "black", size = 10),
    axis.text.x = element_text(color = "black", size = 9, angle = 0, hjust = 1),
    legend.position = "none"
  ) +
  ggtitle(paste0("Female | ", snp_plot, " genotype vs FGFR3 isoform expression"))




p_GTEx_M | p_GTEx_F







################################################################################################
########################################################################## Do lm()
################################################################################################


# Now filter with sex 

sex_to_run <- 2

inputdf <- Mg_data_s %>% filter(SEX == sex_to_run)


# check the sex race age value to set it right 
#inputdf$SEX <- factor(inputdf$SEX, levels = c(1, 2))
inputdf$RACE <- factor(inputdf$RACE)

### Do lm() 

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
    df.tmp <- inputdf[,c(names(inputdf[9:11]), gsub("_add", "", i))] %>%
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
  
  for(j in names(inputdf[9:11])){
  
    
    #j <- "ENST00000481110.7"
    
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
    
    
    
    # MODEL
  
    fmla_un <- as.formula(paste0(j, "~" , i))
    fmla_adj <- as.formula(paste0(j, "~" , i, " + AGE + RACE "))
    
    # log2 TPM + 1
    # fmla_un  <- as.formula(paste0("log2(", j, " + 1) ~ ", i))
    # fmla_adj <- as.formula(paste0("log2(", j, " + 1) ~ ", i, " + AGE + RACE"))
    
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
    #mean_value_1 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 1) %>% pull(j), na.rm = TRUE),2) 
    #mean_value_2 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 2) %>% pull(j), na.rm = TRUE),2) 
    
    
    # create temporary data frame
    df.lm <- data.frame(
      
      snp = i,
      variable = j,
      sex=sex_to_run,
      
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
      
      model = paste0(formula(fmla_un)[2],"~",formula(fmla_un)[3]),
      
      model_adjust = paste0(formula(fmla_adj)[2],"~",formula(fmla_adj)[3]),
      
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

#sum_M <- df.out
#sum_F <- df.out
#Summary_log2TPM_p1 <- rbind(sum_M,sum_F)
#Summary_rawTPM <- rbind(sum_M,sum_F)
#FF <- rbind(Summary_rawTPM,Summary_log2TPM_p1)
# write.table(FF,'/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/plots/GTEx_summ.csv',row.names = F,col.names = T,quote = F,sep = ',')
