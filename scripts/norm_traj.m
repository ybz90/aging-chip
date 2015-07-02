function norm_traj(input_traj)

    % Takes as input a single cell array containing STMs (e.g. traj_export) and outputs a cell array in the same format
    % The output STMSs are normalized as the number of standard deviations above the mean; the distribution of average values for the first cell cycle are taken for all STMs and used to calculate the mean and standard deviations.
    % This normalization process is done for every fluorescent channel in the STMs.

    num_cells = numel(input_traj); %
    sz = size(input_traj{1}); %dimensions of the STMs, as determined by the first one in input_traj
    num_flu = sz(2) - 3; %the number of columns in each STM minus 3 (X; bud frames; metadata) is the # of flu channels

    % Init new traj_export and cell array for storing initial values
    traj_export = cell(1,num_cells);
    init_values = zeros(num_flu,num_cells); %zero matrix of size num_flu rows x num_cells cols

    % Get the average values of each cell's first cell cycle for each channel to populate init_values
    for i = 1:num_cells
        curr_cell = input_traj{i};
        % Get budding times
        budvals = curr_cell(:,end-1);
        X = curr_cell(:,1);
        cycles = []; % get list of budding cell times
        for m = 1:numel(budvals)
            if budvals(m) == 1
                cycles(end+1) = m;
            end
        end
        cycles = sort(cycles);
        for j = 1:num_flu
            curr_flu = curr_cell(1:cycles(1), j+1);
            init_cycle = mean(curr_flu); %get average of first cycle, for current flu channel
            init_values(j,i) = init_cycle; %add average to init_values
        end
    end

    % Calculate mean and std on initial values
    mean_std = zeros(num_flu,2);
    for k = 1:num_flu
        curr_values = init_values(k,:);
        mean_std(k,1) = mean(curr_values);
        mean_std(k,2) = std(curr_values);
    end

    % Normalize each STM according to the mean and standard deviations for each channel
    for l = 1:num_cells
        curr_cell = input_traj{l}; %get old STM
        temp_STM = curr_cell(:,1); %init new, normalized STM; retain first column (time frames)
        for n = 1:num_flu %normalize flu data and add columns to new STM
            old_flu = curr_cell(:,n+1);
            num_data = numel(old_flu);

            curr_mean = mean_std(n,1);
            curr_std = mean_std(n,2);

            temp_flu = zeros(num_data,1);
            for o = 1:num_data %for every fluorescence measurement
                curr_val = old_flu(o);
                norm_val = (curr_val - curr_mean)/curr_std;
                temp_flu(o) = norm_val;
            end

            temp_STM = horzcat(temp_STM,temp_flu); %add normalized column to temp_STM
        end

        % Add bud times and metadata columns to new STM
        temp_STM = horzcat(temp_STM,curr_cell(:,end-1:end));

        % Add temp_STM to traj_export
        traj_export{l} = temp_STM;
    end

    % Export normalized traj_export
    traj_name = input('What do you want to name your output?:','s');
    output_name = [traj_name,'.mat'];
    save(output_name,'traj_export');

end
