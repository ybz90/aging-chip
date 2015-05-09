function plot_data_single(colN,pos_str)

    % plot_data.m takes the trajectories produced by mask_traj.m and combines that with manually curated data for when a cell has died to visualize traces.
    % For now, we are also manually curating the cell cycle count, but once the Whi5 reporter is integrated, this can be automatically determined based on its localization.
    % Plots trajectories for every cell of the input xy positions, along with vertical lines marking the times of cell budding. Requires an input config file containing the start/end frames of each cell's lifespan and the times of budding. This file has filename in the format of xy01_lifespan.txt and format of cell#,start,end,bud1,bud2...

    % Yuan Zhao 05/08/2015


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
        sz = size(all_traj); % get the dimensions of the combined all_traj super array

        % Import manually curated lifespan data for each cell in each position
        % Add lifespan data to all_lifespan 1 x (num pos) cell array, where each lifespan data array has dim:
        % # of cells rows x 3 cols (cell #, lifespan start frame, lifespan end frame)
        lifespan_file = csvread(['xy',pos,'/xy',pos,'_lifespan.txt']);
        all_lifespan{i} = lifespan_file;
    end

    % Figure for containing subplots of traces for each cell
    figure;
    % For every cell in all_traj...
    for j = 1:sz(2)
        [q,r] = quorem(sym(j-1),sym(colN));
        pos_ID = ['xy',pos_str(q+1)]; % quotient of (cell # in all_traj - 1) and (colN per position) + 1 gives index in (positions) of this cell's original xy position ID
        cell_no = r+1; % remainder of (cell # in all_traj - 1) and (colN per position) gives the cell # out of 7 of this cell in its original xy; where r = 0 is cell #1

        % Frames across which the current cell is alive, from manual curation of movies
        lifespan = all_lifespan{q+1}(r+1,3);
        X = all_lifespan{q+1}(r+1,2):lifespan;

        % Initialize current subplot
        gridcol = 4; % four subplots per row
        subplot(ceil(sz(2)/gridcol),gridcol,j); % (rows,cols,current position); rows = numcells/gridcol rounded up
        styles = ['b','r','m','g']; % color and style options for the different channels
        flu_vals = []; % store intensity values for all fluorescent channels

        % Plot the first fluorescent channel
        curr_flu_val = all_traj(X,j,1); % get current cell's fluorescence intensity values across X
        flu_vals{1} = curr_flu_val; % add to flu_vals array
        curr_trace = cell2mat(flu_vals(1));

        ax1 = gca; % get information for first axes, as reference for other flu channel axes
        plot(X,curr_trace,styles(1));
        xlabel(ax1,'Time in frames');
        ylabel(ax1,'GFP');
        cell_title = ['Position ',cell2mat(pos_ID),'; Cell # ',char(cell_no),'; Lifespan: ',num2str(lifespan)];
        title(cell_title);
        hold on

        % Get the cell cycle data from the lifespan file
        % NOTE: Once Whi5 reporter is integrated, cell cycle data can be automated; for now, this too must be added manually
        cycles = all_lifespan{q+1}(r+1,4:end);
        y1 = get(gca,'ylim'); % height of cell cycle bar
        for k = cycles
            line([k k],y1,'Color','k')
        end

        % % Plot every other fluorescent channel...
        % for flu = 2:sz(3)
        %     % Add axis for current flu channel
        %     curr_ax = axes('Position',get(ax1,'Position'), 'Layer','bottom', 'XAxisLocation','top', 'YAxisLocation','right', 'XColor','k', 'YColor',styles(flu));

        %     curr_flu_val = all_traj(X,j,flu);
        %     flu_vals{flu} = curr_flu_val;
        %     curr_trace = cell2mat(flu_vals(flu));

        %     P = plot(X,curr_trace,styles(flu));

        %     axis(curr_ax, 'off', 'tight');
        %     ylabel(curr_ax,'iRFP (nuclear)');
        % end
        % hold off
    end

end