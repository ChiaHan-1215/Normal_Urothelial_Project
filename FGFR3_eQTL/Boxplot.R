# vln box plot for significant gene TPM with 17SNPs from lm result
library(tibble)
library(tidyverse)
library(ggsci)
library(ggrepel)
library(cowplot)
library(ggpubr)
library(rstatix)
library(dplyr)
library(scales)
library(patchwork)

# load TPM and final analysis, using the same number as analysis 
setwd('/Volumes/LTG/Prokunina Lab/Burkitt Lymphoma project/BLGSP analysis/17SNPs eQTL/Chiahan_anaiysis/')

tp <- read.csv('../source files/for analysis/BL_17SNPs_Exp_master_05192024.csv')

# include sample whcih all variable present: Uganda=1, sex, age, EBV 

tp_ug <- tp %>% 
  filter(Country_Uganda_other == 1) %>%
  filter(!is.na(EBV_final)) %>% 
  filter(!is.na(Age_final) & Age_final != "") %>%
  filter(!is.na(sex_final) & sex_final != "") 

#### check correct asign factor or chr or numeric

tp_ug$sex_final <- as.factor(tp_ug$sex_final)
tp_ug$EBV_final <- as.factor(tp_ug$EBV_final)
tp_ug$Country_Uganda_other <-as.factor(tp_ug$Country_Uganda_other)

### Now load the summary table and filter the sig gene and SNPs and percentage of 0 

lout  <- read.csv('17SNPs_lm_TPM_with_TMB_result.csv')
l.sig <-  lout %>% filter(DataSet=='TPM_expression') %>% filter(p.value_adj_sex_age_ebv < 0.057) %>% filter(perc_sample_over_zero > 50)
l.sig$snp <- gsub('_add','',l.sig$snp)

siggene <- l.sig$variable %>% sort() %>% unique()

# try to do plot 
# use LINC01426 as example 
# i <- "LINC01426"

l.plot <- tp_ug %>% select(PID,EBV_final,siggene,grep("rs[0-9]+$",names(tp_ug),value = T)) 
# l.plot <- mt_ug %>% select(PID,EBV_final,SBS5,grep("rs[0-9]+$",names(tp_ug),value = T)) 
# l.plot$SBS5 <- as.numeric(l.plot$SBS5)
#l.plot$rs1888459_md <- car::recode(l.plot$rs1888459," 'AA'='AG' ")
#l.plot$rs1888459_md <- factor(l.plot$rs1888459_md,levels = c("GG","AG"))
#sp <- "rs1888459_md"


# #######################################################################################################################
# #######################################################################################################################
# #######################################################################################################################
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

gene <- siggene[6:9]

plots <- list()
# gene <- "RUNX1"
# sp <- "rs111457485"
# sp <- "rs55799729"
# sp <- "rs1888459"

for (i in gene){

  # i <- "SBS5"
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



# combined_plot <- plots[["RUNX1"]] + plots[["RUNX1_ENST00000455571.5"]] + plots[["RUNX1_ENST00000482318.5"]] + plots[["RPL3P1"]] + plot_layout(ncol = 4)
combined_plot <- plots[["LINC01426"]] + plots[["LINC01426_ENST00000420877.1"]] + plot_layout(ncol = 2)
print(combined_plot)



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

