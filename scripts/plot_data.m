function plot_data(colN,pos_str)

    % run_plot_data.m takes the trajectories produced by mask_traj.m and combines that with manually curated data for when a cell has died to visualize traces.
    % This input config file contains the start/end frames of each cell's lifespan and the times of budding. The filename is of the pattern of xy01_lifespan.txt and has a data format of cell#,start,end,bud1,bud2...

    % Basically, this function is a wrapper for the individual plotting scripts based on it. It horizontally concatenates every position's trajectories into all_traj; it also adds the manually curated lifespan data into a cell matrix, with the data being a cell for every position.

    % Yuan Zhao 05/11/2015


    % Convert input positions (pos from run_analysis.m) to numbers from strings
    positions = [];
    for g = 1:numel(pos_str)
        positions(g) = str2num(pos_str{g});
    end

    % Array for storing all trajectory data across all cells
    all_traj = [];
    % Array for storing manually curated lifespan data
    all_lifespan = cell(1,numel(positions));

    % For every position...
    for i = positions
        % Rename single digit positions to double accordingly
        if numel(i) == 1
            pos = ['0',num2str(i)];
        else
            pos = num2str(i);
        end

        % Horizontally concatenate traj matrices for every position, forming a super array with dimensions:
        % # of frames x # of cells/traps (from all positions) x # of fluorescent channel
        traj_file = ['xy',pos,'/xy',pos,'_traj.mat'];
        load(traj_file);
        all_traj = horzcat(all_traj,traj);

        % Import manually curated lifespan data for each cell in each position
        % Add lifespan data to all_lifespan cell array (1 x num pos), where each position's lifespan data is an array with dim:
        % # of cells rows x 3 cols (cell #, lifespan start frame, lifespan end frame)
        lifespan_file = csvread(['xy',pos,'/xy',pos,'_lifespan.txt']);
        all_lifespan{i} = lifespan_file;
    end

    % Run plotting scripts here

    % Individual subplots for every trajectory
    data_single(colN,pos_str,all_traj,all_lifespan)

    % All trajectories together on a single plot
    %data_all(colN,pos_str,all_traj,all_lifespan)

end
