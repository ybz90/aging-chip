function data_single_meng(pos_str,gridcol,all_traj,all_lifespan,flu_array,label_array)

    % Plots trajectories for every cell of the input xy positions in its own subplot, along with vertical lines marking the times of cell budding.

    % Yuan Zhao 05/18/2015


    % Min/max replicative lifespan of all cells; used to determine a gradient range for plot style color
    max_cycles = cell(1,numel(all_lifespan)); % the highest # of rep life cycles for each position to be plotted; simply take the # of bud frames as the traps with fewer than the max will fill in a 0 anyway
    min_cycles = cell(1,numel(all_lifespan));
    for c = 1:numel(pos_str)
        curr_life = all_lifespan{c};
        max_cycles{c} = numel(curr_life(1,5:end));
    end
    max_cycles = max(cell2mat(max_cycles));

    % Total number of cells to be plotted
    num_cell_all = 0;
    for j = 1:numel(pos_str)
        curr_life = all_lifespan{j};
        sz = size(curr_life);
        num_cell_all = num_cell_all + sz(1);
    end


    % Figure for containing subplots of traces for each cell
    figure;

    %styles = colormap(winter); %define colormap to be used for trajectories to reflect lifespan
    % Red to Black Colormap
    styles = colormap(gray); % first col of gray has 64 rows; black is [0 0 0], red is [1 0 0]
    styles = horzcat(styles(:,1),zeros(64,2)); %add first col of colormap(gray) to 64x2 mat of 0s
    style_num = numel(styles(:,1)); %round(num_cycles/max_cycles)*style_num will translate lifespan of the cell into an index corresponding to the selected colormap

    % Init running total of subplots
    subplot_num = 0;
    
    

    % For every position in pos_str
    for k = 1:numel(pos_str)

        pos = pos_str{k}; % set current position
        curr_traj = all_traj{k}; %load current trajectory file based on position index
        curr_life = all_lifespan{k}; %load lifespan data

        sz = size(curr_life);
        num_cells = sz(1); % number of rows in current position's lifespan data, ie # of cell trajectories to plot
        
        keeptraj=cell(1,num_cells);
        
        % For every cell in the current traj to be plotted
        for l = 1:num_cells

            cell_ID = curr_life(l,1); %the trap in which the current cell is actually located, as not every cell is to be plotted and has info in curr_life
            
            % Initialize current subplot
            subplot(ceil(num_cell_all/gridcol),gridcol,subplot_num+l); % (rows,cols,current position); rows = num_cell_all/gridcol rounded up; current position is the running total of subplots so far, up to this frame, plus the num of the current cell in this frame

            % Frames across which the current cell is alive, from manual curation of movies
            life_end = curr_life(l,4); % third column of appropriate lifespan data
            life_start = curr_life(l,3); % second col
            X = life_start:life_end; % lifespan range

            % Get number of budding cycles from the curr_life data
            %NOTE: Once Whi5 reporter is integrated, cell cycle data can be automated; for now, this too must be added manually
            cycles = curr_life(l,5:end); %budding frames for the current cell
            cycles = cycles(cycles>0); %remove all 0 entries in the cycles array
            num_cycles = numel(cycles); %replicative life span (# of buds) for current cell, after 0s are removed

            % Define current style, based on num budding cycles and colormap from above
            curr_style_idx = ceil((num_cycles/(max_cycles))*style_num);
            curr_style = styles(curr_style_idx,:);


            % Plot the first fluorescent channel in flu_array
            flu_1 = flu_array(1);
            curr_trace = curr_traj(X,cell_ID,flu_1); % get current cell's fluorescence intensity values across X

            keeptraj{l} = [X', curr_trace];
            
            
            % Plot every other fluorescent channel in flu_array...
            for i = 2:numel(flu_array)
                % Add axis for current flu channel

                flu = flu_array(i);
                curr_trace = curr_traj(X,cell_ID,flu);
                
                temptraj = keeptraj{l};
                 temptraj = cat(2,temptraj, curr_trace);
                 keeptraj{l} = temptraj;
                 
%                  tempa = keeptraj_s{l};
%                  tempa = cat(2,tempa, smooth(curr_trace,3));
%                  keeptraj_s{l} = tempa;

%                 P = plot(X,smooth(curr_trace,3),'Color','k','LineWidth',1);
% 
%                 axis(curr_ax, 'off', 'tight');
%                 ylabel(curr_ax,label_array(i));
            end
            
            

            ax1 = gca; % get information for first axes, as reference for other flu channel axes
            
            curr_trace = keeptraj{l}(:,2);
            cell_title = ['xy',pos,', Cell # ',num2str(cell_ID),': Replicative Lifespan = ',num2str(num_cycles)];
            if length(flu_array)==1
                
                plot(X,smooth(curr_trace,3),'Color',curr_style,'LineWidth',1);

                xlabel(ax1,'Time in frames');
                ylabel(ax1,label_array(1));
                title(cell_title);
                hold on

                % Plot cell death type marker; 1 = normal, 2 = late daughter,
                % 3= popped out, 4 = dying round
                if curr_life(l,2) == 1 % square
                    plot(life_end,curr_trace(life_end-life_start+1),'s','MarkerEdgeColor',curr_style,'MarkerFaceColor',curr_style,'Markersize',12,'LineWidth',2)
                end
                if curr_life(l,2) == 2 % triangle
                    plot(life_end,curr_trace(life_end-life_start+1),'>','MarkerEdgeColor',curr_style,'MarkerFaceColor','w','Markersize',10,'LineWidth',2)
                end
                if curr_life(l,2) == 3 % x
                    plot(life_end,curr_trace(life_end-life_start+1),'x','MarkerEdgeColor',curr_style,'MarkerFaceColor','w','Markersize',10,'LineWidth',2)
                end
                if curr_life(l,2) == 4 % o
                    plot(life_end,curr_trace(life_end-life_start+1),'o','MarkerEdgeColor',curr_style,'MarkerFaceColor','w','Markersize',10,'LineWidth',2)
                end

                % Plot cell cycle data as vertical lines
                y1 = get(gca,'ylim'); % height of cell cycle bar
                for k = cycles
                    line([k k],y1,'Color','k','LineStyle','--')
                end


                hold off
            else
                
                [ax, h1, h2]=plotyy(X, keeptraj{l}(:,2),X, keeptraj{l}(:,[3:end]) );
                
                xlabel(ax(1),'Time in frames');
                ylabel(ax(1),label_array(1));
                ylabel(ax(2),label_array(2));
                title(cell_title);
                
                set(h1, 'linewidth',2)
                set(h2, 'linewidth',1.5)
                
                hold on
                axes(ax(1));
                xlim([life_start, life_end]);
                y1 = get(gca,'ylim');
                for k = cycles
                    line([k k],y1,'Color','k','LineStyle','--')
                end
                
                 % Plot cell death type marker; 1 = normal, 2 = late daughter,
                % 3= popped out, 4 = dying round
                if curr_life(l,2) == 1 % square
                    plot(life_end,curr_trace(life_end-life_start+1),'s','MarkerEdgeColor',curr_style,'MarkerFaceColor',curr_style,'Markersize',12,'LineWidth',2)
                end
                if curr_life(l,2) == 2 % triangle
                    plot(life_end,curr_trace(life_end-life_start+1),'>','MarkerEdgeColor',curr_style,'MarkerFaceColor','w','Markersize',10,'LineWidth',2)
                end
                if curr_life(l,2) == 3 % x
                    plot(life_end,curr_trace(life_end-life_start+1),'x','MarkerEdgeColor',curr_style,'MarkerFaceColor','w','Markersize',10,'LineWidth',2)
                end
                if curr_life(l,2) == 4 % o
                    plot(life_end,curr_trace(life_end-life_start+1),'o','MarkerEdgeColor',curr_style,'MarkerFaceColor','w','Markersize',10,'LineWidth',2)
                end
                
                axes(ax(2))
                xlim([life_start, life_end]);
                hold off
            end

        end

        % Increase running total of subplots so far
        subplot_num = subplot_num + num_cells;

    end
    
    

    % % Get all axes handles; set y-limits of all axes
    % axh = findall(gcf,'type','axes');
    % set(axh,'ylim',[0 10000])

end
