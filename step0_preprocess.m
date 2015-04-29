% Create a folder named 'xy position #' and create subfolders for phase,
% threshold, fluorescence, and nuclear channels. Use ImageJ macro to
% pre-process the phase, flu, and nuc channels, and save this into the root
% of the xy position folder. Export these stacks as individual tif images
% into their respective channel/raw subfolders. 

pos = int2str(input('xy position: '));

for channels = {'phase','thresh','flu','nuc'}
    mkdir(strcat('xy',pos,'/',channels{1},'/raw'))
end

% After folder creation, run ImageJ macro for image pre-processing
