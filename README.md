# aging-chip

Pipeline for segmenting, tracking, and analyzing cells in aging chip traps


User Guide:
1. Manually measure the required variables and parameters for ImageJ processing and add them to the run_preprocess_config.txt file (see sample_configs); place this file in working /directory
2. Update the working directory and use the run_preprocess.ijm ImageJ macro to prepare the images for analysis
3. Manually set the lifespan and cell budding times for each position in their respective /directory/xy[POS]_lifespan.txt files (see sample_configs)
4. Use run_analysis.m to generate masks, trajectories, and for visualization analysis


Work in progress:
- Mask generation optimizations (more consistent masks, disappearing nuclei, dilating only if the mother cell is too small)
- Remove small objects below an area representative of cells
- Store non-zero flag for when no cells are detected in the mask
- Declump when the mother is too non-circular, generally reflecting a frame capture during budding
- Use top 20%(?) intensity rather than mean intensity as the measurement of fluorescence
