### Convert files 




```
# In ccad2
ml plink/2.0


plink2 \  
--bfile myFile \  
--recode vcf-iid \  
--keep-allele-order \  
--out myFile


bgzip -c myFile.vcf > myFile.vcf.gz

tabix -p vcf myFile.vcf.gz

```
