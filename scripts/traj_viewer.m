function traj_viewer(traj_export,gridcol,flu_array,label_array)

    % Plots single fluorescent channel trajectory for every cell of the input xy positions in its own subplot, along with vertical lines marking the times of cell budding. Also includes a marker at the end of the trace to signal cell death type.
    % Uses the traj_export output from choose_traj.m, plotting only the selected trajectories.

    % Yuan Zhao 06/22/2015

    % TO DO: ADD CELL DEATH MARKER PLOTTING


    % Total number of cells to be plotted
    num_cell_all = numel(traj_export);


    % Min/max replicative lifespan of all cells; used to determine a gradient range for plot style color
    max_cycles = [];
    for j = 1:num_cell_all
        curr_cell = traj_export{j};
        temp_cycles = curr_cell(:,end-1)>0;
        max_cycles(end+1) = numel(temp_cycles);
    end
    max_cycles = max(max_cycles);

    % Custom colors from rgb database
    green = [0 0.5 0];
    forestgreen = [0.1328 0.5430 0.1328];
    firebrick = [0.6953 0.1328 0.1328];
    crimson = [0.8594 0.0781 0.2344];
    rgbgray = [0.5 0.5 0.5];
    dimgray = [0.4102 0.4102 0.4102];


    % Figure for containing subplots of traces for each cell
    figure;

    %styles = colormap(winter); %define colormap to be used for trajectories to reflect lifespan
    % Red to Black Colormap
    styles = colormap(gray); % first col of gray has 64 rows; black is [0 0 0], red is [1 0 0]
    styles = horzcat(styles(:,1),zeros(64,2)); %add first col of colormap(gray) to 64x2 mat of 0s
    style_num = numel(styles(:,1)); %round(num_cycles/max_cycles)*style_num will translate lifespan of the cell into an index corresponding to the selected colormap


    % For every cell in traj_export...
    for k = 1:num_cell_all

        % Init current subplot
        subplot(ceil(num_cell_all/gridcol),gridcol,k);


        % Get curr_cell data from traj_export
        curr_cell = traj_export{k};

        % First column of curr_cell is the X_range (eg chronological life)
        X = curr_cell(:,1);

        % Get fluorescence data for current cell; for each of the fluorescence channels specified in flu_array, get the appropriate column from 2 to end-2
        flu_plot = cell(1,numel(flu_array));
        for l = 1:numel(flu_array)
            flu = flu_array(l);
            flu_plot{l} = curr_cell(:,flu+1);
        end

        % Get budding times from (end-1)th column; the X-times where values in this column == 1 are budding times
        budvals = curr_cell(:,end-1);
        cycles = [];
        for m = 1:numel(budvals)
            if budvals(m) == 1
                cycles(end+1) = X(m);
            end
        end
        num_cycles = numel(cycles);

        % Title information & death state for subplot
        pos = curr_cell(1,end);
        cell_ID = curr_cell(2,end);


        % Plot subplot
        if pos < 10
            pos = ['0',num2str(pos)];
        else
            pos = num2str(pos);
        end
        cell_title = ['xy',pos,', cell #',num2str(cell_ID),': RLS = ',num2str(num_cycles)];

        % if there is only one fluorescent channel
        if numel(flu_array) == 1
            ax1 = gca;

            % Define current style, based on num budding cycles and colormap from above
            curr_style_idx = ceil((num_cycles/(max_cycles))*style_num);
            curr_style = styles(curr_style_idx,:);

            plot(X,smooth(flu_plot{1},3),'Color',curr_style,'LineWidth',1.5);

            xlabel(ax1,'Time in frames');
            ylabel(ax1,label_array(1));
            title(cell_title);
            hold on

            % Plot cell cycle data as vertical lines
            y1 = get(gca,'ylim'); % height of cell cycle bar
            for k = cycles
                line([k k],y1,'Color','k','LineStyle',':')
            end
        else
            [ax, h1, h2]=plotyy(X, smooth(flu_plot{1},3), X,smooth(flu_plot{2},3) );

            xlabel(ax(1),'Time in frames');
            ylabel(ax(1),label_array(1));
            ylabel(ax(2),label_array(2));
            title(cell_title);

            set(ax,{'ycolor'},{crimson;green})

            set(h1, 'color',crimson, 'linewidth',1.5)
            set(h2, 'color',green, 'linewidth',1.5)

            hold on
            axes(ax(1));
            xlim([X(1), X(end)]);
            y1 = get(gca,'ylim');
            for k = cycles
                line([k k],y1,'Color','k','LineStyle',':')
            end

            axes(ax(2))
            xlim([X(1), X(end)]);
            hold off
        end

    end
end
