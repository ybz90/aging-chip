%% Pipeline wrapper for aging chip code
% Yuan Zhao 05/06/2015

% Input the position number to be registered; to align and account for
% camera movement and drift
pos = int2str(input('Input xy position: '));

% Automatically scan phase folder for the input position to determine the total # of frames
D = dir(['xy',pos,'/c1/']);
imN = length(D(not([D.isdir])));
fprintf('Processing %d frames for position xy %d.\n', imN, str2num(pos));

% (DEBUG) Input the number of frames (times) to be processed
%imN = input('Number of frames to process: ');

%
%% (OPTIONAL) REGISTER ALL FRAMES
% If there is microscape stage shake or drift, register frames
% of all positions to adjust for this movement. Use only if necessary.

register_frames(pos,imN)

%
%% GENERATE MASKS AND CELL TRAJECTORIES

mask_traj(pos,imN,7,2)

%
%% FURTHER ANALYSIS
% put function(s) for plotting and analysis here so the data
