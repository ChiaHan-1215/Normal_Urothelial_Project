### Convert files 

Ref: 

```
plink \  
--bfile myFile \  
--recode vcf-iid \  
--keep-allele-order \  
--out myFile


bgzip -c myFile.vcf > myFile.vcf.gz

tabix -p vcf myFile.vcf.gz

```
