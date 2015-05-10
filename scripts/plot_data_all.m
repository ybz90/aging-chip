function plot_data_single(colN,pos_str)

    %

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

    % Figure containing composite of traces for all cells
    figure;
    % For every cell in all_traj...
    for j = 1:sz(2)
        [q,r] = quorem(sym(j-1),sym(colN));
        pos_ID = ['xy',pos_str(q+1)]; % quotient of (cell # in all_traj - 1) and (colN per position) + 1 gives index in (positions) of this cell's original xy position ID
        cell_no = r+1; % remainder of (cell # in all_traj - 1) and (colN per position) gives the cell # out of 7 of this cell in its original xy; where r = 0 is cell #1

        % Frames across which the current cell is alive, from manual curation of movies
        lifespan = all_lifespan{q+1}(r+1,3);
        X = all_lifespan{q+1}(r+1,2):lifespan;

        styles = ['r','b','m']; % color and style options for the different channels
        flu_vals = []; % store intensity values for all fluorescent channels

        hold on
        % Plot every fluorescent channel...
        for flu = 1%:sz(3)
            curr_flu_val = all_traj(X,j,flu);
            flu_vals{flu} = curr_flu_val;
            curr_trace = cell2mat(flu_vals(flu));

            plot(X,curr_trace,styles(flu))

        end
        hold off

        % Plot fluorescence w/ nuclear marker on the same plot
        %[ax,p1,p2] = plotyy(X,B,X,C,'plot');
        %grid(ax(1),'on')
        %p1.LineWidth = 2;
        %p2.LineWidth = 2;
        %p2.LineSTyle = '--';

        %ylabel(ax(1),'GFP');
        %ylabel(ax(2),'iRFP (nuclear)');
        %xlabel(ax(2),'Time in frames');
        cell_title = ['Position ',cell2mat(pos_ID),'; Cell # ',char(cell_no),'; Lifespan: ',num2str(lifespan)];
        title(cell_title);
    end
end


