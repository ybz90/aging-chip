function plot_data(colN,pos_str)

    % This function is a wrapper for the individual plotting scripts based on it. It loads every position's trajectories into all_traj and the manually curated lifespan data into all_lifespan.

    % Yuan Zhao 05/15/2015


    % Array for storing all trajectory data across all cells
    all_traj = cell(1,numel(pos_str));
    % Array for storing manually curated lifespan data
    all_lifespan = cell(1,numel(pos_str));

    % For every position...
    for i = 1:numel(pos_str)
        pos = pos_str{i};

        % Horizontally concatenate traj matrices for every position, forming a super array with dimensions:
        % # of frames x # of cells/traps (from all positions) x # of fluorescent channel
        traj_file = ['xy',pos,'/xy',pos,'_traj.mat'];
        load(traj_file);
        all_traj{i} = traj;

        % Import manually curated lifespan data for each cell in each position
        % Add lifespan data to all_lifespan cell array (1 x num pos), where each position's lifespan data is an array with dim:
        % # of cells rows x 3 cols (cell #, lifespan start frame, lifespan end frame)
        lifespan_file = csvread(['xy',pos,'/xy',pos,'_lifespan.txt']);
        all_lifespan{i} = lifespan_file;
    end

    % Run plotting scripts here

    % Individual subplots for every trajectory
    data_single(pos_str,8,all_traj,all_lifespan)

    % All trajectories together on a single plot
    %data_all(colN,pos_str,all_traj,all_lifespan)

end
