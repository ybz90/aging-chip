# aging-chip

Pipeline for segmenting, tracking, and analyzing cells in aging chip traps


User Guide:

1. Manually measure the required variables and parameters for ImageJ processing (see sample configs)
2. Update the working directory and use the run_preprocess ImageJ macro to prepare the images for analysis
3. Manually set the lifespan and cell budding times for each position (see sample configs)
4. Use run_analysis to generate masks, trajectories, and visualizations



Work in progress:
- Mask generation optimizations (more consistent masks, disappearing nuclei, etc.)
- Remove small objects below an area representative of cells
- Declump when the mother is too non-circular, generally reflecting a frame capture during budding
- Use top 50%(?) intensity rather than mean intensity as the measurement of fluorescence
