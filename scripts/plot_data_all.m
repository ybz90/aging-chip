function plot_data_all(colN,pos_str,all_traj,all_lifespan)

    % Plots trajectories of every cell of the input xy positions together on a single large plot.

    % Yuan Zhao 05/11/2015


    % Figure containing composite of traces for all cells
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

        flu_vals = []; % store intensity values for all fluorescent channels

        % Define current style, based on lifespan and colormap from above
        curr_style_idx = ceil((lifespan/max_life)*style_num);
        curr_style = styles(curr_style_idx,:);

        hold on
        % Plot every fluorescent channel...
        for flu = 1%:sz(3)
            curr_flu_val = all_traj(X,j,flu);
            flu_vals{flu} = curr_flu_val;
            curr_trace = cell2mat(flu_vals(flu));

            plot(X,curr_trace,'Color',curr_style)

        end
        hold off

        cell_title = 'TITLE';
        title(cell_title);
    end
end


