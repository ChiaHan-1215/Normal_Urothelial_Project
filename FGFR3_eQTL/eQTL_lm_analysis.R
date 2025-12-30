###########################################################
      ## lm analysis ref from BL project #######
###########################################################

library(dplyr)
library(qvalue)
library(lmPerm)

# load dataset

setwd('/Volumes/LTG/Prokunina Lab/Burkitt Lymphoma project/BLGSP analysis/17SNPs eQTL/source files/for analysis/')

## TPM expression of list of genes, (total TPM, iso_tpm)
# need to include a column of how many column is non 0 

tp <- read.csv('BL_17SNPs_Exp_master_05192024.csv')

###### File look like this: The col are gene TPM, and after that is SNP genotype 

#                PID    EBV_final sex_final Age_final Country_Uganda_other   DOP1B  CLIC6   RUNX1    CBR1   MORC3  CHAF1B   SETD4
#1 BLGSP-71-06-00004 ebv-negative         M        11                    1 42.7622 0.3171 12.0381 19.9568 26.4500 33.3083  9.0627
#2 BLGSP-71-06-00095 ebv-negative         M         9                    1  7.2609 2.4477 16.0768  8.5877 12.2717 25.8700  8.9293
#3 BLGSP-71-06-00137 ebv-negative         M        13                    1 22.3748 0.8577 17.5371 16.4104 24.6396 17.4073  9.8266
#4 BLGSP-71-06-00160 ebv-negative         M        19                    1 34.7392 3.0674 10.2850 19.1481 40.9063 19.6882 14.1388
#5 BLGSP-71-08-00021 ebv-negative         M         7                    1 33.8268 0.4979 19.4394 13.5966 57.6424 28.2743 10.6921


tp_ug <- tp %>% 
  filter(Country_Uganda_other == 1) %>%
  filter(!is.na(EBV_final)) %>% 
  filter(!is.na(Age_final) & Age_final != "") %>%
  filter(!is.na(sex_final) & sex_final != "") 


#### check correct asign factor or chr or numeric

tp_ug$sex_final <- as.factor(tp_ug$sex_final)
tp_ug$EBV_final <- as.factor(tp_ug$EBV_final)
tp_ug$Country_Uganda_other <-as.factor(tp_ug$Country_Uganda_other)


# Can make metadata list 

megdata <- list()

megdata[["TPM"]] <- tp_ug
megdata[["TMB_SBS"]] <- mt_ug


####################################################################
####################################################################
################# NOW WE DO lm() ###################################
####################################################################
####################################################################

 
# ####################################################################
         ### TPM/quantile (Z-score) normailzed ####
# ####################################################################

 
inputdf <- megdata$TPM  

# 
df.out <- data.frame()
df_count.summary <- data.frame()
# 




for (i in grep("_add",names(inputdf),value = T)){
  # i <- "rs7279257_add"
  # i <- "rs111457485_add" 
  # i <- "rs16997913_add"
  
  
  function.names <- c("max", "min", "mean", "median", "sd")
  
  for ( k in function.names ) {
    # k <- function.names[1]
    # the all the TPM in all tissue
    # consider dealing with NA
    
    # & inputdf[[i]] != 0
    # names(inputdf[6:135]
    
    df.tmp <- inputdf[,c(names(inputdf[8:16]), gsub("_add", "", i))] %>%
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
    
    # j <- "CBR3.AS1_ENST00000654957.1"
    # j <- names(inputdf[8:16])[1]
    # j <- names(inputdf[6:135])[20]
    
    # 3 model, all sample size should be same, add Permutation column
    
    # unadjust: TPM ~ SNP 
    # adjust: TPM ~ SNP + age + sex
    # adjust for EBV: TPM ~ SNP + age + sex +EBV
    
    # for subsets
    # dynamically generate formula
    fmla_un <- as.formula(paste0(j, "~" , i))
    fmla_adj <- as.formula(paste0(j, "~" , i, " + Age_final + sex_final "))
    fmla_adj_ebv <- as.formula(paste0(j, "~" , i, " + Age_final + sex_final + EBV_final"))
    
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

# write.table(cb,'/Volumes/LTG/Prokunina Lab/Burkitt Lymphoma project/BLGSP analysis/17SNPs eQTL/17SNPs_lm_TPM_with_TMB_result.csv',col.names = T,row.names = F,sep = ',',quote = F)
