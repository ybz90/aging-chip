%% Pipeline wrapper for aging chip code
% Yuan Zhao 05/06/2015

% Add scripts folder to search path
addpath('/Users/yuanz/Git/aging-chip');

% Input the position numbers to be analyzed
pos_input = input('xy positions to analyze {01, 02, ..., 99}: ');

% Convert input positions from ints to strs
pos = {};
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
%% (OPTIONAL) REGISTER ALL FRAMES
% If there is microscape stage shake or drift, register frames of all positions to adjust for this movement. Use only if necessary.

for i = 1:posn
    curr_pos = pos{i};
    fprintf('Registering frames for position xy%d.\n', pos{i});
    register_frames(pos{i},imN)
end

%
%% GENERATE MASKS AND CELL TRAJECTORIES

% Input the number of traps per image
colN = input('Input trap #: ');
% Input the number of fluorescent channels to analyze
fluN = input('Input # of fluorescent channels: ');

for i = 1:posn
    curr_pos = pos{i};
    fprintf('Generating mask and trajectories for position xy%d.\n', pos{i});
    mask_traj(curr_pos,imN,colN,fluN)
end

%
%% VISUALIZATION & ANALYSIS

plot_data_single()

%%

plot_data_all()
