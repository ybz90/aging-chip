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
            all_traj_len(end+1) = traj_end-traj_start; %add to all_traj_len this cell's chronological lifespan
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
    %all_trunc_traj{1}


    % Specify range of window sizes (tau) to test
    tau_range = 1:10;
    tau_eff = [];

    for tau = tau_range
        % Smooth trajectories by moving average of window size tau
        %tau = 1; %debug, no smoothing
        tau_index = 1:shortest_traj+1 - (tau-1); %range will be from start to end minus (tau-1)

        all_tau_traj = cell(1,num_traj); %array for storing trajectories after smoothing

        % for every truncated trajectory
        for m = 1:num_traj
            temp_traj = []; %build smoothed trajectory in this matrix
            prev_traj = all_trunc_traj{m}; %unsmoothed trajectory

            for n = tau_index %in the smoothed index range
                frames_to_avg = []; %frames within the window for every index in tau_index
                for o = n:n + tau-1
                    %[n,o]
                    frames_to_avg(end+1) = prev_traj(o);
                end
                curr_avg = mean(frames_to_avg);
                temp_traj(end+1) = curr_avg;

            end

            all_tau_traj{m} = temp_traj;
        end
        %all_tau_traj{1}
        %all_tau_traj{3} == (all_trunc_traj{3})' %debug; check equivalence with tau = 1


        % Adjust curves via average of all traj last window
        all_last = [];
        for p = 1:num_traj %for every smoothed trajectory
            curr_traj = all_tau_traj{p};
            all_last(end+1) = curr_traj(end); %take its final window value
        end
        %calculate average of last window for all trajectories; add or subtract the entire trajectories for each one to normalize their end window value to this mean
        avg_last = mean(all_last)
        %std_last = std(all_last)

        all_adj_traj = cell(1,num_traj); % store trajectories after subtraction of delta
        for q = 1:num_traj;
            curr_traj = all_tau_traj{q}; % for every smoothed traj
            delta = curr_traj(end) - avg_last; % find the delta from its end window and the overall mean for end windows
            all_adj_traj{q} = curr_traj - delta; % subtract this delta from all values in array
        end


        % Emsemble average of end-adjusted trajectory start windows
        % This section takes the ensemble average of the tau start points of each trajectory (moving window of size tau) and finds the ensemble standard deviation as well.
        all_start = [];
        for r = 1:num_traj %for every smoothed trajectory
            curr_traj = all_adj_traj{r}; % of the trajectories adjusted for end window value sync
            all_start(end+1) = curr_traj(1); %take its first window value
        end
        avg_start = mean(all_start);
        std_start = std(all_start)


        % Unsilencing event detection
        % Based on the threshold, representing a multiplier of the above ensemble standard deviation, if the final window average of the trajectory is greater than the threshold + that trajectory's start window average, increase one_count by 1. If it does not meet detection criterion, increase zero_count by 1.
        % Use this to find the detection efficiency of the reporter for the current window size tau and threshold.
        threshold = 2;

        one_count = 0;
        zero_count = 0;

        for s = 1:num_traj
            curr_traj = all_adj_traj{s};
            curr_start = curr_traj(1); %start window value of current traj
            target = curr_start + threshold*std_start; %threshold value = curr_start + x*std; end window must exceed this to count as 1
            if avg_last > target
                one_count = one_count + 1;
            else
                zero_count = zero_count + 1;
            end
        end

        % calculate efficiency
        efficiency = one_count/(one_count+zero_count)

        % store in tau_eff matrix
        tau_eff(end+1) = efficiency;

    end

    tau_range
    tau_eff

end
