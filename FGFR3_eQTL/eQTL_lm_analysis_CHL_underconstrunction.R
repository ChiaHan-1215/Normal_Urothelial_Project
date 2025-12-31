library(dplyr)
library(qvalue)
library(lmPerm)



df_fgfr <- read.csv('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Parse_scRNA-seq/Sample_Genotyping/FGFR3_qn.snp.Test.csv')

df_fgfr <- df_fgfr %>% 
  filter(!is.na(Predicted_Sex) & Predicted_Sex != "") %>%
  filter(!is.na(Ancestry) & Ancestry != "") 


#### check correct asign factor or chr or numeric

df_fgfr$Predicted_Sex <- as.factor(df_fgfr$Predicted_Sex)
df_fgfr$Ancestry <-as.factor(df_fgfr$Ancestry)

#### Make SNP as 0,1,2 

# Iterate through SNP columns
for (i in names(df_fgfr)[5:26]) {
  
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


df_fgfr <- df_fgfr[,c(names(df_fgfr)[1:4],setdiff(names(df_fgfr),names(df_fgfr)[1:4]) %>% sort())]

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




for (i in grep("rs",names(inputdf),value = T)){
  # i <- "rs1374468"
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
  
  for(j in names(inputdf[6:135])){
    
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
    
    
    # make 0 as 1, not sure we need here
    
    filtered_0_1_2_only <- filtered_0_1_2_only %>% mutate(!!i := ifelse(.[[i]] == 0, 1, .[[i]]))
    
    # fit lm model
    fit_un <- lm(fmla_un, filtered_0_1_2_only)
    fit_adj <- lm(fmla_adj, filtered_0_1_2_only)
    fit_adj_ebv <- lm(fmla_adj_ebv, filtered_0_1_2_only)
    
    
    ########################################
    # run permutation #
    ########################################
    
    prm <- lmp(fmla_un, data = filtered_0_1_2_only, perm = "Prob")
    prm_adj <- lmp(fmla_adj, data = filtered_0_1_2_only, perm = "Prob")
    prm_adj_ebv <- lmp(fmla_adj_ebv, data = filtered_0_1_2_only, perm = "Prob")
    
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
      
      mean_value_0_and_1=mean_value_1,
      mean_value_2=mean_value_2,
      
      
      
      # add lm() result
      p.value = coef(summary(fit_un))[2,4],
      p.prm = coef(summary(prm))[2,3],
      beta = coef(summary(fit_un))[2,1],
      StdEr = round(coef(summary(fit_un))[2,2], digits = 4),
      
      p.value_adj_sex_age = coef(summary(fit_adj))[2,4],
      p.prm_adj_sex_age = coef(summary(prm_adj))[2,3],
      beta_adj_sex_age = coef(summary(fit_adj))[2,1],
      StdEr_adj_sex_age = round(coef(summary(fit_adj))[2,2], digits = 4),
      
      p.value_adj_sex_age_ebv = coef(summary(fit_adj_ebv))[2,4],
      p.prm_adj_sex_age_ebv = coef(summary(prm_adj_ebv))[2,3],
      beta_adj_sex_age_ebv =  coef(summary(fit_adj_ebv))[2,1],
      StdEr_adj_sex_age_ebv = round(coef(summary(fit_adj_ebv))[2,2], digits = 4)
      
    )
    
    
    # bind rows of temporary data frame to the results data frame
    
    df.out <- rbind(df.out, df.lm)
    
    
  }
}



################################################################################################################################################
################################################################################################################################################


# now cb the dataset 

# f1 <- read.csv('~/Desktop/BL_project/17SNPs_lm_TPM_result.csv')
# f2 <- read.csv('~/Desktop/BL_project/17SNPs_lm_TPM_summary.csv')

set1 <- df.out
set2 <- df.out.tmb


set1 <- set1 %>%
  mutate(DataSet = "TPM_expression") %>%
  select(1:which(names(set1) == "variable"), DataSet, (which(names(set1) == "variable") + 1):ncol(set1))


set2 <- set2 %>%
  mutate(DataSet = "Mutation") %>%
  select(1:which(names(set2) == "variable"), DataSet, (which(names(set2) == "variable") + 1):ncol(set2))


cb <- rbind(set1,set2)


#  cb.summ <- inner_join(f2,df_count.tmb.summary)
