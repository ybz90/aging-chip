%% Pipeline wrapper for aging chip code
% Yuan Zhao 06/23/2015

%% Initialization

% This code section adds aging-chip scripts to path and specifies positions to be analyzed.

% Add scripts folder to search path
addpath('/Users/yuanz/Git/aging-chip/scripts');
addpath('/Users/yuanz/Git/aging-chip/dev'); %folder for scripts in development

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

% This sections runs the mask_traj code on each position specified above individually and generates masks and trajectories. It will run in parallel.

tic
% Input the number of fluorescent channels to analyze
fluN = input('Input # of fluorescent channels: ');

parfor i = 1:posn
    curr_pos = pos{i};
    fprintf('Generating mask and trajectories for position xy%d.\n', str2double(curr_pos));
    gen_mask_traj(curr_pos,imN,fluN)
end
toc

%
%% Trajectory Processing

% The following code processes the output trajectories from above, based on the positions input above and xy##_lifespan.txt files specified for each positions.
% It loads every position's trajectories into all_traj and the manually curated lifespan data into all_lifespan.

% Array for storing all trajectory data across all cells
all_traj = cell(1,numel(pos));
% Array for storing manually curated lifespan data
all_lifespan = cell(1,numel(pos));

% For every position...
for i = 1:numel(pos)
    curr_pos = pos{i};

    % Horizontally concatenate traj matrices for every position, forming a super array with dimensions:
    % # of frames x # of cells/traps (from all positions) x # of fluorescent channel
    traj_file = ['xy',curr_pos,'/xy',curr_pos,'_traj.mat'];
    load(traj_file);
    all_traj{i} = traj;

    % Import manually curated lifespan data for each cell in each position
    % Add lifespan data to all_lifespan cell array (1 x num pos), where each position's lifespan data is an array with dim:
    % # of cells rows x 3 cols (cell #, lifespan start frame, lifespan end frame)
    lifespan_file = csvread(['xy',curr_pos,'/xy',curr_pos,'_lifespan.txt']);
    all_lifespan{i} = lifespan_file;
end

disp('all_traj and all_lifespan loaded for input positions.');


% Choose the trajectories to be exported
% Input fluorescent channels to plot (c2 = 1, c3 = 2, etc.)
flu_array = input('Fluorescent channels to plot [1 2 ...]: ');

% Input fluorescent channel labels
%label_array = input('Fluorescent channel labels {'GFP','irFP',etc.}');
label_array = {'mCherry','GFP'};

gridcol = 5; %indicate number of columns of subplots
choose_traj(pos,gridcol,all_traj,all_lifespan,flu_array,label_array)

%
%% Run Analysis & Visualization Code

% Plot trajectories
%traj_viewer(traj_export,5,flu_array,label_array)

% Trajectory analysis scripts
traj_stat(norm_export,1,norm_export,2,1)

%
