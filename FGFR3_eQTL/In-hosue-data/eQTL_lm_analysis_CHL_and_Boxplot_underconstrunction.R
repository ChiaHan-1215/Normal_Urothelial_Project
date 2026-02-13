library(dplyr)
library(tximport)
library(BSgenome.Hsapiens.UCSC.hg38)
library(EnsDb.Hsapiens.v86)
library(xlsx)
library(edgeR)





####################################################################
####################################################################
################# NOW WE DO lm() ###################################
####################################################################
####################################################################


# ####################################################################
### TPM/quantile (Z-score) normailzed ####
# ####################################################################


inputdf <- megdata$FGFR3  
# inputdf <- Indf
# 
inputdf <- read.csv('../FGFR3_isoform_TMM_INT.snp_500k.csv')
str(inputdf)

inputdf$Predicted_Sex <- factor(inputdf$Predicted_Sex)
inputdf$Ancestry <- factor(inputdf$Ancestry)

df.out <- data.frame()
df_count.summary <- data.frame()
# 




for (i in grep("_add",names(inputdf),value = T)){
  # i <- "rs62286992_add"
  
  
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
  
  for(j in names(inputdf)[c(10,16,11)]){
    
    # j <- "ENST00000481110_FGFR3IIIc"
    
    
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
    
    
    # dynamically generate formula
    fmla_un <- as.formula(paste0(j, "~" , i))
    fmla_adj <- as.formula(paste0(j, "~" , i, " + Predicted_Sex +  Ancestry"))
    fmla_adj_2 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex +  Ancestry +",i,paste0("*Predicted_Sex")))
    fmla_adj_3 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex + RIN_score + AFR + EUR + ASN"))
    fmla_adj_4 <- as.formula(paste0(j, "~" , i, " + Predicted_Sex + RIN_score + AFR + EUR + ASN + ",i,paste0("*Predicted_Sex")))
    
    
    ######################################################
    ######### make 0 as 1, not sure we need here #########
    # filtered_0_1_2_only <- filtered_0_1_2_only %>% mutate(!!i := ifelse(.[[i]] == 0, 1, .[[i]]))
    ######################################################
    
    # fit lm model
    fit_un <- lm(fmla_un, filtered_0_1_2_only)
    fit_adj <- lm(fmla_adj, filtered_0_1_2_only)
    fit_adj_2 <- lm(fmla_adj_2, filtered_0_1_2_only)
    fit_adj_3 <- lm(fmla_adj_3, filtered_0_1_2_only)
    fit_adj_4 <- lm(fmla_adj_4, filtered_0_1_2_only)
    
    ########################################
    # run permutation #
    ########################################
    
    #prm <- lmp(fmla_un, data = filtered_0_1_2_only, perm = "Prob")
    #prm_adj <- lmp(fmla_adj, data = filtered_0_1_2_only, perm = "Prob")
    
    # Calculate the number of non-zero samples in j where i is not NA
    
    # use this since the lm() is done based on the sample inculde GT
    non_zero_count <- sum(filtered_0_1_2_only[[j]] != 0 & !is.na(filtered_0_1_2_only[[j]]))
    
    # non_zero_count <- sum(inputdf[[j]] != 0 & !is.na(inputdf[[j]]))
    
    # Calculate mean values for SNP dosages 1 and 2
    # mean_value_1 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 1) %>% pull(j), na.rm = TRUE),2) 
    # mean_value_2 <- round(mean(filtered_0_1_2_only %>% filter(.[[i]] == 2) %>% pull(j), na.rm = TRUE),2) 
    
    
    
    
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
      
      #p.prm = tryCatch({coef(summary(prm))[2,3]}, error=function(e){
      #  print(paste0('error of NA'))
      #  return(NA)}),
      
      # Main effects (wrapped in tryCatch to prevent loop breaks)
      p.value = tryCatch(coef(summary(fit_un))[2,4], error = function(e) NA),
      beta = tryCatch(round(coef(summary(fit_un))[2,1], 4), error = function(e) NA),
    
      #StdEr = round(coef(summary(fit_un))[2,2], digits = 4),
      
      p.value_adj = coef(summary(fit_adj))[2,4],
      beta_adj = coef(summary(fit_adj))[2,1],
      # StdEr_adj_sex_age = round(coef(summary(fit_adj))[2,2], digits = 4),
      # p.prm_adj_sex_age = coef(summary(prm_adj))[2,3],
      p.value_adj_2_int = tryCatch(coef(summary(fit_adj_2))[grep(paste0("^",i, ":"), rownames(coef(summary(fit_adj_2)))), "Pr(>|t|)"],error = function(e) NA)[1],
      
      p.value_adj_3 = tryCatch(coef(summary(fit_adj_3))[2,4], error = function(e) NA),
      beta_adj_3 = tryCatch(coef(summary(fit_adj_3))[2,1], error = function(e) NA),
      
      p.value_adj_4_int = tryCatch(coef(summary(fit_adj_4))[grep(paste0("^",i, ":"), rownames(coef(summary(fit_adj_4)))), "Pr(>|t|)"],error = function(e) NA)[1]
      
    )
    
    
    # bind rows of temporary data frame to the results data frame
    
    df.out <- rbind(df.out, df.lm)
    
    
  }
}



##################################
# Test the Boxplot 
##################################

# 3 SNPs p-value is sig, try to do plot 
# rs1374468_add
# rs2854915_add
# rs11930034_add


# Use the data used for lm() analysis 
# inputdf

# vln box plot for significant gene TPM with 17SNPs from lm result
library(ggplot2)
library(ggsci)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(rstatix)
library(dplyr)
library(scales)
library(patchwork)

# try to do plot 
# use LINC01426 as example 
# i <- "LINC01426"

sig_snp <- c('rs1374468','rs2854915','rs11930034')

l.plot <- inputdf %>% dplyr::select(Sample_Name,FGFR3,Ancestry,rs1374468,rs1063743_add,rs2854915,rs2854915_add,rs11930034,rs11930034_add) 

# ############################ FOR changing GT order ####################################################################
# 
# 
# for (i in 1:nrow(l.sig)){
# 
# # i <- 11
# 
# sp <- l.sig$snp[i]
# gene <- l.sig$variable[i]
# 
# l.plot_order <- l.plot %>%
#   filter(!is.na(l.plot[[sp]]) & l.plot[[sp]] != "")
# 
# s1 <- c(l.sig$geno_1[i],l.sig$geno_2[i])
# 
# if(sp == "rs7279257"){
#   l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = c('TT','AT','AA'))
#   } else {l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = s1)}
# 
# 
# 
# p <-l.plot_order %>%
#   ggplot(aes(x = .data[[sp]], y = .data[[gene]], fill = .data[[sp]])) +
#   geom_violin(width = 1, lwd = 0.2, colour = "black",
#               show.legend = NA, inherit.aes = TRUE) +
#   # geom_violin(trim = FALSE) +
#   scale_fill_manual(values = c("#f4a582", "#92c5de", "#0571b0"),drop=F) + ## modified
#   
#   geom_boxplot(
#     width=0.15,lwd = 0.5,fill="white", position = position_dodge(width = 1),
#     outlier.shape=NA) + 
#   
#   geom_point(shape = 1, size = 0.75, colour = "black",
#              position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) +
#   # Add the mean point
#   stat_summary(fun = mean, geom = "point",
#                position = position_dodge(width = 1),
#                shape = 19, size = 1.5, colour = "red") + theme_classic() + 
#   theme(axis.title.x=element_blank(),
#         axis.text.y = element_text(color = "black",size = 10),
#         axis.text.x = element_text(color = "black",size = 10), 
#         axis.title.y =element_blank(),legend.position = "none") + 
#   ggtitle(paste(sp,'with expression of',gene,'\n p-value:',round(l.sig$p.value_adj_sex_age_ebv[i],3),', beta:',round(l.sig$beta_adj_sex_age_ebv[i],3))) 
# 
# 
# 
# print(p)
# 
# }
# 
# 
# 

########### plot rs11145 and rs224  ########### 

# cahnge order of rs111457485
sig_snp
plots <- list()
gene <- "FGFR3"
sp <- sig_snp[3]

for (i in gene){
  
  # i <- "FGFR3"
  counts <- l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "") %>%  group_by(.data[[sp]]) %>% summarise(n=n())
  x_labels <- paste(counts[[sp]], "\n(n=", counts$n, ")", sep = "")
  # y_max <- max(l.plot[[i]], na.rm = TRUE) + 0.25
  
  
  p <-l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "") %>% 
    ggplot(aes(x = .data[[sp]], y = .data[[i]], fill = .data[[sp]])) +
    geom_boxplot(
      width=0.5,lwd = 0.5, 
      outlier.shape=NA) + 
    
    scale_fill_manual(values = c("#F38491", "#6ECFF6", "#4197EC"),drop=F) +
    
    geom_point(shape = 1, size = 1.5, colour = "black",
               position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1)) +
    # Add the mean point
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 1),
                 shape = 19, size = 1.5, colour = "red") + theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color = "black",size = 10),
          axis.text.x = element_text(color = "black",size = 8,angle = 45, hjust = 1), 
          axis.title.y =element_blank(),legend.position = "none",
          plot.title = element_text(size = 10, face = "plain", color = "black")) + scale_x_discrete(labels = x_labels) + #ylim(0, 45) +
    ggtitle(paste0("Expression of ",i," vs ",sp))
  
  
  #print(p)
  plots[[i]] <- p
  
}

################################################################################
# CUT OFF
#################################################################################



######################################## for rs224, need to change the GT 
################################################################################
################################################################################

sp <- "rs2242780"
sp <- "rs7279257_md"
gene <- siggene[2:3]
# i <- gene[1]
plots <- list()

l.plot_order <- l.plot %>% filter(!is.na(.data[[sp]]) & .data[[sp]] != "")
l.plot_order$rs7279257_md <- car::recode(l.plot_order$rs7279257, " 'TT'='AT'   ")
#l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = c("GG","CG"),ordered = T)
l.plot_order[[sp]] <- factor(l.plot_order[[sp]],levels = c("CC","CT"),ordered = T)

for (i in gene){
  counts <- l.plot_order %>%  group_by(.data[[sp]]) %>% summarise(n=n())
  x_labels <- paste(counts[[sp]], "\n(n=", counts$n, ")", sep = "")
  # y_max <- max(l.plot[[i]], na.rm = TRUE) + 0.25
  
  
  p <-l.plot_order %>% 
    ggplot(aes(x = .data[[sp]], y = .data[[i]], fill = .data[[sp]])) +
    geom_boxplot(
      width=0.5,lwd = 0.5, 
      outlier.shape=NA) + 
    
    scale_fill_manual(values = c("#F38491", "#6ECFF6", "#4197EC"),drop=F) +
    
    geom_point(shape = 1, size = 1.5, colour = "black",
               position = position_jitterdodge(jitter.width = 0.2,dodge.width = 1)) +
    # Add the mean point
    stat_summary(fun = mean, geom = "point",
                 position = position_dodge(width = 1),
                 shape = 19, size = 1.5, colour = "red") + theme_classic() + 
    theme(axis.title.x=element_blank(),
          axis.text.y = element_text(color = "black",size = 10),
          axis.text.x = element_text(color = "black",size = 8,angle = 45, hjust = 1), 
          axis.title.y =element_blank(),legend.position = "none",
          plot.title = element_text(size = 10, face = "plain", color = "black")) + scale_x_discrete(labels = x_labels) +
    ggtitle(paste0("Expression of ",i," vs ",sp))
  
  
  #print(p)
  #plots[[i]] <- p
  plots[[paste0(i,"__rs5")]] <- p
}

combined_plot <- plots[["LINC01426"]] + plots[["LINC01426_ENST00000420877.1"]] +  plots[["LINC01426__rs5"]] + plots[["LINC01426_ENST00000420877.1__rs5"]] + plot_layout(ncol = 4)
print(combined_plot)
