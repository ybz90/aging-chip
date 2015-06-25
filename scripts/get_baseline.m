function get_baseline(traj_export,flu_ch)

    % Takes as input a traj_export set of cell trajectories and exports the 'baseline', an average of each trajectory's initial value for a given fluorescence channel.
    % To be used on traj_export outputs from choose_traj.m.

    % Yuan Zhao 06/23/2015


    % Total number of cells to be plotted
    num_cell_all = numel(traj_export);

    % Initial values for each trajectory in traj_export
    initial = [];

    % For every cell in traj_export...
    for k = 1:num_cell_all

        % Get curr_cell data from traj_export
        curr_cell = traj_export{k};

        % Initial value of curr_cell for flu_ch
        i_0 = curr_cell(1,flu_ch+1);

        initial(end+1) = i_0;

    end

    mean(initial)

end
