function combine_same_traj(traj_array,flu_keep)

	% Inputs are a set of cell arrays in the format of traj/norm_export and an array indicating which flu channel to keep in that particular array
	% then we will take those cell arrays, discard all unwanted flu channels, and combine them all into a super array
	% the intended use is to combine the same type of data (e.g. all NTS data) from multiple experiments into one array, provided conditions are similar or the data has been normalized with get_baseline and then norm_traj

	% yz 6/26

	traj_all = {};

	% For each set of traj to be 'cleaned'
	for i = 1:numel(traj_array)

		curr_traj = traj_array{i}; % current traj set, of the form of traj/norm_export
		curr_flu = flu_keep{i}; % channel to keep in this traj_set

		num_cell = numel(curr_traj);

		cleaned = cell(1,num_cell); %temp to store curr_traj after retaining only the desired flu ch

		for j = 1:num_cell

			% current cell; ie matrix within curr_traj cell array
			curr_cell = curr_traj{j};

			% keep first col of matrix (ie X), desired flu ch col, and last two col of matrix (ie bud frames, metadata) for current cell
			new_cell = horzcat(curr_cell(:,1),curr_cell(:,curr_flu+1),curr_cell(:,end-1:end));

			cleaned{j} = new_cell;

		end

		traj_all = horzcat(traj_all,cleaned);

	end

	traj_name = input('What do you want to name your output?:','s');
	output_name = [traj_name,'.mat'];
	save(output_name,'traj_all');

end
