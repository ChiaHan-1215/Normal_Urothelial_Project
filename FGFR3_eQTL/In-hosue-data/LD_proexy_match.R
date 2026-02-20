# Use LDlinkR to extract SNP proxy

# Date: 02092026

library(LDlinkR)


# EUR hg38 coord 

my_proxies <- LDproxy(snp = "rs2896518", 
                      pop = "EUR", 
                      r2d = "r2", 
                      token = "c22a071b5bda",genome_build = "grch38_high_coverage",win_size = 500000
)


my_proxies$Coord <- gsub(':','_',my_proxies$Coord)

# Download LD proxy output, extract the GT from In-house data and see how many GT are exist in the porxy output


library(dplyr)

setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/LD_proxy_files/')

# read in house data
Indf <- read.csv('../FGFR3_isoform_TMM_INT.snp_500k_v2.csv')

Our_target <- names(Indf) %>% grep('chr4_[0-9+]+$|rs[0-9]+$',.,value = T) %>% as.data.frame()
names(Our_target) <- "Coord"

target_for_rsid <- Our_target %>% filter(grepl("^rs[0-9]+$", Coord))
names(target_for_rsid) <- "RS_Number"


cb_for_pos <- inner_join(Our_target %>% filter(!grepl("^rs[0-9]+$", Coord)),my_proxies,by = 'Coord')
cb_for_pos <- cb_for_pos[,c(2,1,3,4:11)]
cb_for_rsid <- inner_join(target_for_rsid,my_proxies,by = 'RS_Number')

ff <- rbind(cb_for_pos,cb_for_rsid)


# ###########################################
# ############ CUT OFF OLD ##################
# ###########################################
# 
# 
# 
# 
# 
# # Download LD proxy output, extract the GT from In-house data and see how many GT are exist in the porxy output
# # change the LD proxy output location from hg19 to hg38
# # Date: 02092026
# 
# library(dplyr)
# 
# setwd('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/LD_proxy_files/')
# # proxy <- read.table('proxy_rs2896518_ALL.pop_R2_No_collapse_FORGEdb_500000bp.txt',header = T,stringsAsFactors = F)
# 
# 
# # exract and save the rsID, and use UCSC table browser to input the SNP IDs and outptu the location 
# # https://genome.ucsc.edu/cgi-bin/hgTables
# # Parameter Setting: Assembly: hg38, Grooup: Variation, Track: dbSNP 155, Table: All dbSNP(155)(dbSnp 155)
# # and upload the rsid list in the "Define region of interest --> UPload list"
# # chose as csv to save
# 
# #write.table(proxy %>% filter(RS_Number != '.') %>% select(RS_Number),'proxy_rs2896518_SNPlist.txt',col.names = F,row.names = F,sep = '\t',quote = F)
# 
# 
# # read the lifted file
# USCS_lifted <- read.csv('UCSC_liftover_hg38.csv')
# USCS_lifted <- USCS_lifted %>% filter(X.chrom == 'chr4')
# names(USCS_lifted)[1:4] <- c('chrom','hg38_start','hg38_end',"RS_Number")
# 
# USCS_lifted <- USCS_lifted[,c(4,1,2,3,5,6)]
# 
# cb <- left_join(proxy,USCS_lifted,by='RS_Number')
# names(cb)[2] <- "hg19_Coord"
# 
# 
# cb$SNP_hg38_pos <- paste0(cb$chrom,"_",cb$hg38_end)
# 
# # read in house data
# Indf <- read.csv('../FGFR3_isoform_TMM_INT_with_SEX_RACE_for_lm.snp.csv')
# 
# Our_target <- names(Indf) %>% grep('chr4_[0-9+]+$|rs[0-9]+$',.,value = T) %>% as.data.frame()
# names(Our_target) <- "SNP_hg38_pos"
# 
# target_for_rsid <- Our_target %>% filter(grepl("^rs[0-9]+$", SNP_hg38_pos))
# names(target_for_rsid) <- "RS_Number"
# 
# cb_for_pos <- inner_join(Our_target %>% filter(!grepl("^rs[0-9]+$", SNP_hg38_pos)),cb,by = 'SNP_hg38_pos')
# cb_for_pos <- cb_for_pos[,c(2,1,3,4:12)]
# 
# cb_for_rsid <- inner_join(target_for_rsid,cb,by = 'RS_Number')
# cb_for_rsid <- cb_for_rsid[,c(1,17,2:11)]
# 
# ff <- rbind(cb_for_pos,cb_for_rsid)



# Extract the lm result to extract the SNPs that is in proxy data

res_lm <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/Victor_Normal_Urothelial_project/Project_FGFGR3/lm_FGFR3_iso_TMM_INT_vs_SNPs_adj_int_v2.tsv')
names(res_lm)[1] <- "RS_Number"
res_lm_rsONLY <- res_lm[(grepl('rs',res_lm$RS_Number)),]
res_lm_POS <- res_lm[!(grepl('rs',res_lm$RS_Number)),]

COMBINE_rs <- inner_join(ff,res_lm_rsONLY,by = 'RS_Number')
# 
# SET1 <- ff %>% filter(!(RS_Number %in% unique(COMBINE_rs$RS_Number)))
# names(SET1) 
# SET1$RS_Number <- SET1$Coord
# 
# COMBINE_Pos <- inner_join(SET1,res_lm_POS,by = "RS_Number")
#Final_data <- rbind(COMBINE_rs,COMBINE_Pos)

# Now seperate risk and protective allele
#  rs2896518-A risk
library(stringr)

COMBINE_rs <- COMBINE_rs %>%
  mutate(
    Correlated_risk    = str_extract(Correlated_Alleles, "(?<=A=)[^,]+"),
    Correlated_protect = str_extract(Correlated_Alleles, "(?<=G=).+")
  )  %>% relocate(Correlated_risk,Correlated_protect,.before = variable)


