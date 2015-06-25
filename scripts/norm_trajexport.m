function norm_trajexport(traj_export,flu_ch,scale)

    % Takes as input a traj_export, flu_ch, and the calculated scaling factor and multiples the appropriate flu_ch by the scaling factor.

    % Yuan Zhao 06/23/2015


    % Total number of cells to be plotted
    num_cell_all = numel(traj_export);

    % Normalized cell array to replace traj_export
    norm_export = cell(1,num_cell_all);

    % For every cell in traj_export...
    for k = 1:num_cell_all

        % Get curr_cell data from traj_export
        curr_cell = traj_export{k};

        % Break current cell's traj array into slices (1) preceding flu_ch column to scale; (2) flu_ch column to scale; and (3) following flu_ch column to scale
        pre = curr_cell(:,1:flu_ch);
        flu = curr_cell(:,flu_ch+1);
        post = curr_cell(:,flu_ch+2:end);

        % scale flu channel data
        flu = flu*scale;

        % Concatenate array together again
        scaled_cell = horzcat(pre,flu,post);

        % Add scaled cell to norm_export
        norm_export{k} = scaled_cell;

    end

traj_name = input('What do you want to name your output?: ','s');
output_name = [traj_name,'.mat'];
save(output_name,'norm_export');

end
