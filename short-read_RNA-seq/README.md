#### The RNA-seq notes of urothelial project
- Goal
  
  The short read RAN-seq results of normal urothelial project

  Furture directions:
    - run GSEA analysis/GO term to determine cell biological bacground
  


  
------------------
- Data
  
  fastq files (read length: 150bp) delivederd from CGR with QC sheets etc. the location is in T-drive:
  
    - `T:\DCEG\Projects\DataDelivery\Prokunina\TP0325-RS7-Urothelial-Samples-RNA-seq`
  
  sample: 116 sample described in `Sample_table.md`. as for the RNA-seq QC sheet, `114` sample.

  When checking, the missing 3 sample: `SD347567,SD347590, SD347628` , and `SC917163` is not include in original sample sheet. what's that?

  BAM files
  
  

  
  
------------------
- Code

  following GTEx RNA-seq piepline: https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq
  or searching CGR github to see if there are good built piepline

  Following CGR pipeline, do STAR indexing and alignment, the script name `star_index_and_align.md`
  
  STAR alignmetn
  RSEM to generated TPM
  RseQC to QC
  
- Results
-------------------
