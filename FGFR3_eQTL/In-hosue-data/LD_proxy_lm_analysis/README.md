
## Goal:

start from LD proxy result, find which allele is corrected with risk/protect, then check allele is flip or not.

### The new approach:

	- LD=proxy result => select GT from it => then see mismatch or not => flip => set risk/pro allele
	- => use risk as 2 and portect as 0 => and then do lm() 


#### For isoform TPMs 

```R

## Goal:

start from LD proxy result, find which allele is corrected with risk/protect, then check allele is flip or not.

### The new approach:

	- LD=proxy result => select GT from it => then see mismatch or not => flip => set risk/pro allele
	- => use risk as 2 and portect as 0 => and then do lm() 

### Do boxplot for R2 over 0.5
