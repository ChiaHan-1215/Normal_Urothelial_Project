### Convert files 




```
# In ccad2

# or plink 1.9?
# ml plink/1.9


ml plink/2.0
ml bgzip/1.15 tabix/1.15 



plink2 \  
--bfile myFile \  
--recode vcf-iid \  
--keep-allele-order \  
--out myFile


bgzip -c myFile.vcf > myFile.vcf.gz

tabix -p vcf myFile.vcf.gz

```
