function plot_data_single(colN,pos_str,all_traj,all_lifespan)

    % For now, we are also manually curating the cell cycle count, but once the Whi5 reporter is integrated, this can be automatically determined based on its localization.
    % Plots trajectories for every cell of the input xy positions in its own subplot, along with vertical lines marking the times of cell budding.

    % Yuan Zhao 05/11/2015


    % Figure for containing subplots of traces for each cell
    figure;

    % Min/max lifespan of all cells; used to determine a gradient range for plot style color
    max_life = [];
    min_life = [];
    for c = 1:numel(all_lifespan)
        max_life = vertcat(max_life,all_lifespan{c}(:,3));
        %min_life = vertcat(min_life,all_lifespan{c}(:,2));
    end
    max_life = max(max_life);
    %styles = colormap(winter); %define colormap to be used for trajectories to reflect lifespan
    % Red to Black Colormap
    styles = colormap(gray); % first col of gray has 64 rows; black is [0 0 0], red is [1 0 0]
    styles = horzcat(styles(:,1),zeros(64,2)); %add first col of colormap(gray) to 64x2 mat of 0s
    style_num = numel(styles(:,1)); %round(lifespan/max_life)*style_num will translate lifespan of the cell into an index corresponding to the selected colormap

    % For every cell in all_traj...
    sz = size(all_traj); % get the dimensions of the combined all_traj super array
    for j = 1:sz(2)
        [q,r] = quorem(sym(j-1),sym(colN));
        pos_ID = ['xy',pos_str(q+1)]; % quotient of (cell # in all_traj - 1) and (colN per position) + 1 gives index in (positions) of this cell's original xy position ID
        cell_no = r+1; % remainder of (cell # in all_traj - 1) and (colN per position) gives the cell # out of 7 of this cell in its original xy; where r = 0 is cell #1

        % Frames across which the current cell is alive, from manual curation of movies
        life_end = all_lifespan{q+1}(r+1,3); % third column of appropriate lifespan data, based on pos ID
        life_start = all_lifespan{q+1}(r+1,2); % second col
        X = life_start:life_end; % lifespan range
        lifespan = life_end-life_start;

        % Initialize current subplot
        gridcol = 5; % four subplots per row
        subplot(ceil(sz(2)/gridcol),gridcol,j); % (rows,cols,current position); rows = numcells/gridcol rounded up
        flu_vals = [];  % store intensity values for all fluorescent channels

        % Define current style, based on lifespan and colormap from above
        curr_style_idx = ceil((lifespan/(max_life))*style_num);
        curr_style = styles(curr_style_idx,:);

        % Plot the first fluorescent channel
        curr_flu_val = all_traj(X,j,1); % get current cell's fluorescence intensity values across X
        flu_vals{1} = curr_flu_val; % add to flu_vals array
        curr_trace = cell2mat(flu_vals(1));

        ax1 = gca; % get information for first axes, as reference for other flu channel axes
        plot(X,curr_trace,'Color',curr_style);
        xlabel(ax1,'Time in frames');
        ylabel(ax1,'GFP');
        cell_title = ['Position ',cell2mat(pos_ID),'; Cell # ',char(cell_no),'; Lifespan: ',num2str(lifespan)];
        title(cell_title);
        hold on

        % Plot cell cycle data as vertical lines; get the cell cycle data from the lifespan file
        % NOTE: Once Whi5 reporter is integrated, cell cycle data can be automated; for now, this too must be added manually
        cycles = all_lifespan{q+1}(r+1,4:end);
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

end
