# aging-chip

Pipeline for segmenting, tracking, and analyzing cells in aging chip traps

step0_init_preprocess.ijm:
Takes xy01c1.tif > xy01_c1_t.tif, xy01_c1_thr_t.tif, etc.

step1_maskgen_register.m:
Takes xy01_c1_t.tif > xy01_c1_reg_t.tif, etc. 

step2_maskgen_mask.m:
Takes xy01_c1_thr_reg_t.tif, xy01_c3_reg_t.tif > xy01_mask_raw_t.tif

step3_maskgen_clean.m:
Takes xy01_mask_raw_t.tif, xy01_c1_reg_t.tif, xy01_c1_thr_reg_t.tif, xy01_c2_reg_t.tif > xy01_mask_clean_t.tif, xy01_mask_combined_t.tif