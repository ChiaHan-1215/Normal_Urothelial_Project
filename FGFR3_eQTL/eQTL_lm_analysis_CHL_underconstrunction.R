
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
