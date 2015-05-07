# aging-chip

Pipeline for segmenting, tracking, and analyzing cells in aging chip traps


Work in progress:
- Mask generation optimizations (more consistent masks, disappearing nuclei, dilating only if the mother cell is too small)
- Remove small objects below an area representative of cells
- Store non-zero flag for when no cells are detected in the mask
- Declump when the mother is too non-circular, generally reflecting a frame capture during budding
- Use top 20%(?) intensity rather than mean intensity as the measurement of fluorescence
