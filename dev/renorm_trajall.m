function renorm_trajall(input_traj,div)

    % rescales a cell array of STMs with a single flu channel
    % temp code thats not really meant for keeping, just to rescale 06/26 data from baseline of 1000 to 200

    %div = 5; %divide by 5; we're rescaling 1000 to 200

    num_cells = numel(input_traj);
    traj_all = cell(1,num_cells); %output traj_all, same var name and size as input

    for i = 1:num_cells
        curr_cell = input_traj{i};

        temp = curr_cell(:,1); %get first column

        flu = curr_cell(:,2)/div; %rescale flu data

        temp = horzcat(temp,flu); %and rescaled flu data
        temp = horzcat(temp,curr_cell(:,end-1:end)); % add last two columns

        traj_all{i} = temp;

    end

    traj_name = input('What do you want to name your output?:','s');
    output_name = [traj_name,'.mat'];
    save(output_name,'traj_all');

end
