% NO LONGER NECESSARY; FUNCTIONALITY IMPLEMENTED IN IMAGEJ MACRO

% Create a folder named 'xy position #' and create subfolders for phase,
% threshold, fluorescence, and nuclear channels. Use ImageJ macro to
% pre-process the phase, flu, and nuc channels, and save this into the root
% of the xy position folder. Export these stacks as individual tif images
% into their respective channel/raw subfolders. 


% Input xy position #s of images to be stitched
pos_input = input('xy positions to initialize folders for {01, 02, ..., 99}: ');

% Convert input positions from ints to strs
pos = {};
for k = 1:length(pos_input) % for each position input
    pos{k} = int2str(pos_input{k}); % add its xy pos name to pos {}
    if length(pos{k}) == 1 % for single digit positions, append 0 in front to make it two digits, ie 1 -> 01
        pos{k} = ['0',pos{k}];
    end
end
posn = length(pos);

% create folders for each xy_pos, and subfolders for ph, thr, flu, nuc, and
% /raw subfolder in each
for i = 1:posn
    j = pos{i};
    for channels = {'phase','thresh','flu','nuc'}
        mkdir(['xy',j,'/',channels{1},'/raw'])
    end
end


% After folder creation, run ImageJ macro for image pre-processing