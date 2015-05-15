function data_single(colN,pos_str,all_traj,all_lifespan)

    % For now, we are also manually curating the cell cycle count, but once the Whi5 reporter is integrated, this can be automatically determined based on its localization.
    % Plots trajectories for every cell of the input xy positions in its own subplot, along with vertical lines marking the times of cell budding.

    % Yuan Zhao 05/15/2015

    % TO DO: Set consistent axes for comparison
    % TO DO: Instead of plotting all cells 1-7, only plot the ones in the lifespan.txt config file


    % Figure for containing subplots of traces for each cell
    figure;

    % Min/max replicative lifespan of all cells; used to determine a gradient range for plot style color
    max_cycles = cell(1,numel(all_lifespan)); % the highest # of rep life cycles for each position to be plotted; simply take the # of bud frames as the traps with fewer than the max will fill in a 0 anyway
    min_cycles = cell(1,numel(all_lifespan));
    for c = 1:numel(pos_str)
        curr_life = all_lifespan{c};
        max_cycles{c} = numel(curr_life(1,4:end));
    end
    max_cycles = max(cell2mat(max_cycles));

    %styles = colormap(winter); %define colormap to be used for trajectories to reflect lifespan
    % Red to Black Colormap
    styles = colormap(gray); % first col of gray has 64 rows; black is [0 0 0], red is [1 0 0]
    styles = horzcat(styles(:,1),zeros(64,2)); %add first col of colormap(gray) to 64x2 mat of 0s
    style_num = numel(styles(:,1)); %round(num_cycles/max_cycles)*style_num will translate lifespan of the cell into an index corresponding to the selected colormap


    % Total number of cells to be plotted
    num_cell_all = 0;
    for j = 1:numel(pos_str)
        curr_life = all_lifespan{j};
        sz = size(curr_life);
        num_cell_all = num_cell_all + sz(1);
    end

    % Init running total of subplots
    subplot_num = 0;

    % For every position in pos_str
    for k = 1:numel(pos_str)

        pos = pos_str{k}; % set current position
        curr_traj = all_traj{k}; %load current trajectory file based on position index
        curr_life = all_lifespan{k}; %load lifespan data

        sz = size(curr_life);
        num_cells = sz(1); % number of rows in current position's lifespan data, ie # of cell trajectories to plot

        % For every cell in the current traj to be plotted
        for l = 1:num_cells

            cell_ID = curr_life(l,1); %the trap in which the current cell is actually located, as not every cell is to be plotted and has info in curr_life

            % Initialize current subplot
            gridcol = 8; % four subplots per row
            subplot(ceil(num_cell_all/gridcol),gridcol,subplot_num+l); % (rows,cols,current position); rows = num_cell_all/gridcol rounded up; current position is the running total of subplots so far, up to this frame, plus the num of the current cell in this frame
            flu_vals = [];  % store intensity values for all fluorescent channels

            % Frames across which the current cell is alive, from manual curation of movies
            life_end = curr_life(l,3); % third column of appropriate lifespan data
            life_start = curr_life(l,2); % second col
            X = life_start:life_end; % lifespan range

            % Get number of budding cycles from the curr_life data
            %NOTE: Once Whi5 reporter is integrated, cell cycle data can be automated; for now, this too must be added manually
            cycles = curr_life(l,4:end); %budding frames for the current cell
            cycles = cycles(cycles>0); %remove all 0 entries in the cycles array
            num_cycles = numel(cycles); %replicative life span (# of buds) for current cell, after 0s are removed

            % Define current style, based on num budding cycles and colormap from above
            curr_style_idx = ceil((num_cycles/(max_cycles))*style_num)
            curr_style = styles(curr_style_idx,:);


            % Plot the first fluorescent channel
            curr_flu_val = curr_traj(X,cell_ID,1); % get current cell's fluorescence intensity values across X
            flu_vals{1} = curr_flu_val; % add to flu_vals array
            curr_trace = cell2mat(flu_vals(1));

            ax1 = gca; % get information for first axes, as reference for other flu channel axes
            plot(X,curr_trace,'Color',curr_style,'LineWidth',1);
            xlabel(ax1,'Time in frames');
            ylabel(ax1,'GFP');
            cell_title = ['Position ',pos,', Cell # ',num2str(cell_ID),'; Replicative Lifespan: ',num2str(num_cycles)];
            title(cell_title);
            hold on

            %Plot cell cycle data as vertical lines
            y1 = get(gca,'ylim'); % height of cell cycle bar
            for k = cycles
                line([k k],y1,'Color','k','LineStyle','--')
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

        % Increase running total of subplots so far
        subplot_num = subplot_num + num_cells;

    end

    % % Get all axes handles; set y-limits of all axes
    % axh = findall(gcf,'type','axes');
    % set(axh,'ylim',[0 10000])

end
