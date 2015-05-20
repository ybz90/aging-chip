%% Pipeline wrapper for aging chip code
% Yuan Zhao 05/13/2015

%% Initialization

% Add scripts folder to search path
addpath('/Users/yuanz/Git/aging-chip/scripts');

% Input the position numbers to be analyzed
pos_input = input('xy positions to analyze {01, 02, ..., 99}: ');

% Convert input positions from ints to strs
pos = cell(1,numel(pos_input));
for k = 1:length(pos_input) % for each position input
    pos{k} = int2str(pos_input{k}); % add its xy pos name to pos {}
    if length(pos{k}) == 1 % for single digit positions, append 0 in front to make it two digits, ie 1 -> 01
        pos{k} = ['0',pos{k}];
    end
end
posn = length(pos);

% Automatically scan the first position's phase folder to determine the total # of frames output by the ImageJ macro
D = dir(['xy',pos{1},'/c1/']);
imN = length(D(not([D.isdir])));
fprintf('Processing %d frames.\n', imN);

% (DEBUG) Input the number of frames (times) to be processed
%imN = input('Number of frames to process: ');

%
%% Generate masks and cell trajectories

tic 

% % Input the number of traps per image
% colN = input('Input # of traps: ');
colN = 7;

% Input the number of fluorescent channels to analyze
fluN = input('Input # of fluorescent channels: ');

parfor i = 1:posn
    curr_pos = pos{i};
    fprintf('Generating mask and trajectories for position xy%d.\n', str2num(curr_pos));
    mask_traj(curr_pos,imN,colN,fluN)
end

toc
%
%% Visualization

% Input fluorescent channels to plot (c2 = 1, c3 = 2, etc.)
flu_array = input('Fluorescent channels to plot [1 2 ...]: ');

% Input fluorescent channel labels
%label_array = input('Fluorescent channel labels {GFP,irFP,etc.}'); 
label_array = {'GFP','mCherry'};

plot_data(colN,pos,flu_array,label_array)

%