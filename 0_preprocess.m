% Create a folder named 'xy position #' and create subfolders for phase,
% threshold, fluorescence, and nuclear channels. Use ImageJ macro to
% pre-process the phase, flu, and nuc channels, and save this into the root
% of the xy position folder. Export these stacks as individual tif images
% into their respective channel/raw subfolders. 


%
% Pre-processing steps; use ImageJ macro (choose an appropriate range of frames; 
% based on phase, change angle of all images to horizontal; subtract background 
% from all flu channels; create threshold images from phase; remove cells/features
% above and below channels of interest in threshold images; export 3 channel 
% stacks and 1 threshold stack)

% By convention: c1 is phase, c2 is fluorescence, c3 is nuclear marker
% c1_pha is the threshold image based on phase
% After registering, the files are renamed as eg. xy02c1_pha_A_txxxx.tif,
% xy02c1_A_txxxx.tif
%


pos = int2str(input('xy position: '));

for channels = {'phase','thresh','flu','nuc'}
    mkdir(strcat('xy',pos,'/',channels{1},'/raw'))
end