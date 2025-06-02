#########################################################
# geting GTEX ID and target SNPs from subsetted vcf file#
#########################################################

library(vcfR)
library(tidyr)
library(reshape2)
library(stringr) 
library(textclean)
library(xlsx)
# snp.list <- data.frame(rsid = c("TERT_TR2",
#                                 "haplotype_2SNP",
#                                 "haplotype_3SNP",
#                                 "rs56345976",
#                                 "rs33961405",
#                                 "rs10069690",
#                                 "rs2242652",
#                                 "rs2853669",
#                                 "rs7712562",
#                                 "rs2735940",
#                                 "rs2736100",
#                                 "rs2853677",
#                                 "rs7705526"),
#                        SNP.pick_hg38 = c("chr5:1275400:G:A",
#                                          "chr5:1276752:G:A",
#                                          "chr5:1276753:G:A",
#                                          "chr5:1276758:G:A",
#                                          "chr5:1277462:G:A",
#                                          "chr5:1279675:C:T",
#                                          "chr5:1279913:G:A",
#                                          "chr5:1295234:A:G",
#                                          "chr5:1295957:A:G",
#                                          "chr5:1296371:A:G",
#                                          "chr5:1286401:C:A",
#                                          "chr5:1287079:G:A",
#                                          "chr5:1285859:C:A" ))

spl <- readxl::read_xlsx('~/Desktop/chr15_CHRNA5_bladder_project/CHRNA5_markers_new.xlsx')

vcf.file.tmp <- read.vcfR("/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/CHRNA5_project/GTEx_v9_WGS_PSMA4_C5.vcf")

# subset VCF for T2/haplotypes
vcf.file.set <- vcf.file.tmp[ which(getPOS(vcf.file.tmp) %in% gsub("15:", "", spl$hg38)), ]

# get and format genotypes, set "return.alleles = TRUE" can switch 0|1 to A|G
vcf_genotypes.tmp <- extract.gt(
  vcf.file.set,
  element = "GT",
  mask = FALSE,
  as.numeric = FALSE,
  return.alleles = TRUE,
  IDtoRowNames = TRUE,
  extract = TRUE,
  convertNA = TRUE
)

# transform dataframe row/column
geno_df.tmp <- as.data.frame(t(vcf_genotypes.tmp))
# use match function, but %in% also work well.
# %in% give T/F logical output while match gives vector position 
# names(geno_df.tmp) <- snp.list$rsid[ match(gsub("^chr5_|_[A-Z]_.+", "", names(geno_df.tmp)), gsub("chr5:|:[A-Z]", "", snp.list$SNP.pick_hg38)) ]
geno_df.tmp$ID_person <- row.names(geno_df.tmp)
geno_df.tmp <- geno_df.tmp[,c(9,1:8)]

# remove "/" or "|"  between genotype
geno_df.tmp <- data.frame(lapply(geno_df.tmp, function(x){ gsub("\\/", " ", x)}))
# write.table(geno_df.tmp,'~/Desktop/chr15_CHRNA5_bladder_project/newset.SNP.GT.for.GTEx.C5andP4.csv',row.names = F,col.names = T,quote = F,sep = ',')

# to remove the last -0003 etc from GTEx_ID
# FOR EXAMPLE: 
#ids <- c("GTEX-1117F-0003", "GTEX-111VG-0004", "GTEX-ABCDE-1234")
# remove the last “-XXXX” chunk (hyphen plus any non-hyphen chars to end)
# geno_df.tmp$GTEx_ID <- sub("-[^-]+$", "", geno_df.tmp$ID_person)

print(donor_ids)
#> [1] "GTEX-1117F"  "GTEX-111VG"  "GTEX-ABCDE"


