function traj_stat(pos_str,all_traj,all_lifespan,tau,threshold)

    % DESCRIPTION
    % tau = window size
    % threshold = # x standard deviation
    flu = 1;

    % Yuan Zhao 06/04/2015


    % Initialization
    % This section will process all of the trajectories of interest by truncating them to the same length (the length of the shortest valid trajectory), all ending at their death points. In other words, each trajectory will be shortened to the index range of (death time - shortest traj len):(death time).
    % Additionally, we 'normalize' all trajectories to their final death point value by averaging them all (this death 'point' is actually an ensemble average of the final points of each trajectory, defined by window size tau), and shifting the entire trajectory up or down to reach this average.

    % Determine the number of trajectories (sample size),
    all_traj_len = []; %chronological lifespans

    for c = 1:numel(pos_str) %for every position
        curr_life = all_lifespan{c}; %load the xy##_lifespan.txt for this position
        for d = 1:numel(curr_life(1:end,1)) %for every row (cell) in this position's lifespan config
            traj_end = curr_life(d,4);
            traj_start = curr_life(d,3);
            all_traj_len = horzcat(all_traj_len,[traj_end-traj_start]); %add to all_traj_len this cell's chronological lifespan
        end
    end

    %all_traj_len

    num_traj = numel(all_traj_len);
    shortest_traj = min(all_traj_len);


    % Use shortest_traj length to truncate all trajectories, ending at the traj_end/death time
    all_trunc_traj = cell(1,num_traj);

    count = 1;
    for k = 1:numel(pos_str)
        pos = pos_str{k}; % set current position
        curr_traj = all_traj{k}; %load current trajectory file based on position index
        curr_life = all_lifespan{k}; %load xy##_lifespan.txt lifespan data for this position

        sz = size(curr_life);
        num_cells = sz(1); % number of rows in current position's lifespan data, ie # of cell trajectories to plot

        for l = 1:num_cells
            cell_ID = curr_life(l,1); % current cell's ID
            traj_end = curr_life(l,4); % end point of current cell's chronological lifespan

            curr_trunc_ind = traj_end - shortest_traj:traj_end; %indices of truncated trajectory
            a = curr_traj(curr_trunc_ind,cell_ID,flu); %slice truncated trajectory from curr_traj based on curr_trunc_ind
            all_trunc_traj{count} = a;

            count = count + 1;
        end
    end

    %all_trunc_traj


    % Emsemble average
    % This section takes the ensemble average of the tau start points of each trajectory (moving window of size tau) and finds the ensemble standard deviation as well.


    % Unsilencing event detection
    % Based on the threshold, representing a multiplier of the above ensemble standard deviation, if the final window average of the trajectory is greater than the threshold + that trajectory's start window average, increase one_count by 1. If it does not meet detection criterion, increase zero_count by 1.
    % Use this to find the detection efficiency of the reporter for the current window size tau and threshold.



    % PLOT?

    % input tau array? run above code for all tau in tau array? plot v tau? for tau>1, try multi start points ~ error bar? for every threshold, we have a different efficiency vs tau curve


end
