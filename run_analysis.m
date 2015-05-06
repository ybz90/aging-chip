%% INPUT CONFIGURATION

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
%% GENERATE MASKS FROM NUCLEAR MARKER

mask_gen(pos,imN)

%
%% TRAJECTORIES
% NOTE: VERIFY MENG'S CODE FOR NOW, GET SOME TRAJECTORIES. THINK OF BETTER WAYS TO IMPLEMENT THIS CODE. ALWAYS EXPORT THIS DATA .MAT BUT KEEP IN MEM TOO FOR LATER ANALYSIS? MAYBE JUST REIMPORT TO BE SAFE IF MULTIPEL ANALYSIS SCRIPTS ARE RUN

autotrack_columns_trajdata_simple_v1(imN,pos,7,2)

%
%% FURTHER ANALYSIS
% put function(s) for plotting and analysis here so the data
