function plot_data_all(colN,pos_str)

    %

    % Yuan Zhao 05/10/2015


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
        % Add lifespan data to all_lifespan cell array (1 x num pos), where each position's lifespan data is an array with dim:
        % # of cells rows x 3 cols (cell #, lifespan start frame, lifespan end frame)
        lifespan_file = csvread(['xy',pos,'/xy',pos,'_lifespan.txt']);
        all_lifespan{i} = lifespan_file;
    end

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
    styles = colormap(winter); %define colormap to be used for trajectories to reflect lifespan
    style_num = numel(styles(:,1)); %round(lifespan/max_life)*style_num will translate lifespan of the cell into an index corresponding to the selected colormap

    % For every cell in all_traj...
    for j = 1:sz(2)
        [q,r] = quorem(sym(j-1),sym(colN));
        pos_ID = ['xy',pos_str(q+1)]; % quotient of (cell # in all_traj - 1) and (colN per position) + 1 gives index in (positions) of this cell's original xy position ID
        cell_no = r+1; % remainder of (cell # in all_traj - 1) and (colN per position) gives the cell # out of 7 of this cell in its original xy; where r = 0 is cell #1

        % Frames across which the current cell is alive, from manual curation of movies
        lifespan = all_lifespan{q+1}(r+1,3);
        X = all_lifespan{q+1}(r+1,2):lifespan;

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

        cell_title = ['Position ',cell2mat(pos_ID),'; Cell # ',char(cell_no),'; Lifespan: ',num2str(lifespan)];
        title(cell_title);
    end
end


