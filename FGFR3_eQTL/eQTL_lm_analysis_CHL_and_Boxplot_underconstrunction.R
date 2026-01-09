
### Load VcfR 
library(vcfR)
library(biomaRt)

# ==============================================================================
# 1) LOAD VCF AND EXTRACT DATA
# ==============================================================================

# load cov

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
geno_df.tmp <- geno_df.tmp[,c("ID_person",setdiff(names(geno_df.tmp),"ID_person"))]
geno_df.tmp <- data.frame(lapply(geno_df.tmp, function(x){ gsub("\\/", " ", x)}))
geno_df.tmp$GTEx_ID <- sub("-[^-]+$", "", geno_df.tmp$ID_person)



# load the sample isofom count


iso <- read.delim('/Volumes/ifs/DCEG/Branches/LTG/Prokunina/GTEx_data/project_FGFR3/FGFR3_iso.txt')
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

library(edgeR)

dge <- DGEList(counts = mani_bladder_with_iso)

# calculate normalization factors
dge <- calcNormFactors(dge, method = "TMM")

# get normalized counts
normalized_counts <- cpm(dge)
head(normalized_counts[1:10,1:10])
dim(normalized_counts)

#

inv_norm_transform <- function(x) {
  qnorm((rank(x, na.last = "keep") - 0.5) / sum(!is.na(x)))
}

# apply inverse normal transform to each gene
final_normalized <- t(apply(normalized_counts, 1, inv_norm_transform))
final_normalized <- as.data.frame(final_normalized)
head(final_normalized[1:5,1:10])

