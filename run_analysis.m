%% STEP 0: CONFIG
% Set folder path to the working directory; add this script's folder to the search path.

% Input the position number to be registered; to align and account for
% camera movement and drift
pos = int2str(input('Input xy position: '));

D = dir(['xy',pos,'/c1/']);
imN = length(D(not([D.isdir])));
fprintf('Processing %d frames for position xy %d.\n', imN, str2num(pos));

%imN = input('Number of frames to process: ');

%
%% STEP 1: (OPTIONAL) REGISTER ALL FRAMES
% If there is microscape stage shake or drift, register frames
% of all positions to adjust for this movement. Use only if necessary.

register_frames(pos,imN)

%
%% STEP 2: GENERATE MASKS FROM NUCLEAR MARKER
% NOTE: GOAL IS TO MERGE THESE TWO MASK GEN STEPS; PERHAPS EVENTUALLY, EVEN COMBINE WITH TRAJECTORY s.t. AS SOON AS MASK IS GENERATED, THOSE OBJECTS ARE ANALYZED AND TRACKS; LEAVE AN EXPORT OPTION AS A DEBUG CHOICE THAT IS GENERALLY COMMENTED OUT

mask_gen(pos,imN)

%
%% TRAJECTORIES
% NOTE: VERIFY MENG'S CODE FOR NOW, GET SOME TRAJECTORIES. THINK OF BETTER WAYS TO IMPLEMENT THIS CODE. ALWAYS EXPORT THIS DATA .MAT BUT KEEP IN MEM TOO FOR LATER ANALYSIS? MAYBE JUST REIMPORT TO BE SAFE IF MULTIPEL ANALYSIS SCRIPTS ARE RUN

autotrack_columns_trajdata_simple_v1(imN,pos,1,2)

%
%% FURTHER ANALYSIS
% put function(s) for analysis here so the data
