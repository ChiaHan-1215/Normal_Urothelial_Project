#### The RNA-seq notes of urothelial project
- Goal
  
  The short read RAN-seq results of normal urothelial project 
- Data
  
  fastq files delivederd from CGR with QC sheets etc. the location is in T-drive:
  
  `T:\DCEG\Projects\DataDelivery\Prokunina\TP0325-RS7-Urothelial-Samples-RNA-seq`
  
  sample: 116 sample described in `Sample_table.md`.
  
  The RNA-seq QC sheet, `114` sample.
  
  read length: 150bp

- Code
  
  following GTEx RNA-seq piepline: https://github.com/broadinstitute/gtex-pipeline/tree/master/rnaseq
  or searching CGR github to see if there are good built piepline

  Following CGR pipeline, do STAR indexing and alignment, the script name `star_index_and_align.md`
  
  STAR alignmetn
  RSEM to generated TPM
  RseQC to QC
  
- Results
