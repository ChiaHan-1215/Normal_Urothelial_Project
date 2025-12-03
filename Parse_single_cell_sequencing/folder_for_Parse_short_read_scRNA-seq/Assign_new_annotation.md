## Goal: Assign annotitaion from UCSC cell browser: https://cells.ucsc.edu/?ds=tabula-sapiens+by-organ+bladder&meta=free_annotation

- add this annotation to our rds file
- also add the compartment annotation

## file location and updates:

file location: `/data/leec20/parse_single_cell/2024_output/MEGA_cb_all_75_bladder_sample/all-sample/Per_sample_filtering_and_cellAnnted/`

- The original loading file: `BP_FINAL_MERGEed_1M_parse_withcellanno_Intered.rds`
- BPcell location relocated to new folder, need to relink
- add the SCTtransform to the rds file
- add new annotation


## Code: 

Code are located in the Biowulf: `/data/leec20/project_FGFR3_GADD4_parse_single_cell/`


## Output files

- The file wuth SCT: `BP_FINAL_MERGEed_1M_parse_withcellanno_SCT.rds` 
- The file with SCT and corrected BPcell location: `BP_FINAL_MERGEed_1M_parse_withcellanno_SCT_corrected_BP.rds`
- The file with SCT, corrected BPcell location, and annotation: `BP_FINAL_MERGEed_1M_parse_withcellanno_SCT_corrected_BP_NewAnno.rds`
