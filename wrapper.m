% Input the position number to be registered; to align and account for
% camera movement and drift
pos = int2str(input('Input xy position: '));

% (DEBUG) Input the number of frames (times) to be processed
imN = input('Number of frames to process: '); 

%% STEP 1: REGISTER ALL FRAMES (OPTIONAL)
% (Optional) If there is microscape stage shake or drift, register frames
% of all positions to adjust for this movement. 

register_frames(pos,imN)
 
%% STEP 2: GENERATE MASKS FROM NUCLEAR MARKER

mask_gen(pos,imN)

%% STEP 3: CLEAN MASK

mask_clean(pos,imN)

%% TRAJECTORIES

autotrack_columns_trajdata_simple_v1(imN,pos,1,2)