function combine_traj(traj_array,flu_keep)

	% Inputs are two arrays; one corresponding to a set of cell arrays that contain matrices in standard trajectory matrix (STM) format, and the other is a list of the respective fluorescent channels to be kept in each of the input cell arrays
	% For each pair of STM-containing cell array and respective fluorescent channel, this code removes all other fluorescent channels from the matrices in that cell array
	% Then, all of these 'cleaned' cell arrays are horizontally concatenated to combine all of the STMs into a single, large cell array
	% The intend use is to combine compatible data from multiple experiments (e.g. all NTS2 from identical conditions in multiple experiments; or NTS2 from various conditions after 'pseudo-normalization') into a single array of STMs with just one fluorescent dat a channel

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
