function plot_data

    % plot_data.m takes the trajectories produced by mask_traj.m and combines that with manually curated data for when a cell has died to visualize traces. For now, we are also manually curating the cell cycle count, but once the Whi5 reporter is integrated, this can be automatically determined based on its localization.

    % Yuan Zhao 05/07/2015


    % READ FROM EXT COFNIG FILE S.T. FOR EVERY CELL WE KNOW WHAT ITS LIFE
    % IS FRAME RANGE AS MANUALLY COUNTED
    % WHEN THAT CELL IS BEING PLOTTED, IT'S DATA WILL BE 1:LIFE instead of
    % 1:sz(1), which is the full frame #


    colN = 7; % there are seven traps / cells per position
    positions = 27:29; %read this from configfile


    % Array for storing all trajectory data across all cells
    all_traj = [];

    % Horizontally concatenate traj matrices for every position, forming a super array with dimensions:
    % # frames x # cells/traps (from all positions) x # fluorescent channel
    for i = positions
        pos = num2str(i);
        traj_file = ['xy',pos,'/xy',pos,'_traj.mat']
        load(traj_file);
        all_traj = horzcat(all_traj,traj);
        sz = size(all_traj); % get the dimensions of the combined all_traj super array
    end


    % Figure for containing subplots of traces for each cell
    figure;
    X_all = 1:sz(1); % define global axes

    % For every cell in all_traj...
    for j = 1:sz(2)
        flu_vals = []; % store intensity values for all fluorescent channels

        [q,r] = quorem(sym(j-1),sym(colN));
        pos_ID = ['xy',num2str(positions(q+1))]; % quotient of (cell # in all_traj - 1) and (colN per position) + 1 gives index in (positions) of this cell's original xy position ID
        cell_no = r+1; % remainder of (cell # in all_traj - 1) and (colN per position) gives the cell # out of 7 of this cell in its original xy; where r = 0 is cell #1

        % Frames across which the current cell is alive, from manual curation of movies
        lifespan = 400;
        X = 1:lifespan; % SET TO EXACT LIFE BASED ON MANUAL CHECKING; WE CAN MAP TO THAT # BY LOOKING AT POSITION XY AND 1-7 CELL/COLN;

        % For every fluorescent channel...
        for flu = 1:sz(3)
            curr_flu_val = all_traj(X,j,flu); % get current cell's fluorescence intensity values across X
            flu_vals{flu} = curr_flu_val; % add to flu_vals array
        end

        C = cell2mat(flu_vals(2));
        B = cell2mat(flu_vals(1));

        gridcol = 4; % four subplots per row
        ax = subplot(ceil(sz(2)/gridcol),gridcol,j); % (rows,cols,current position); rows = numcells/gridcol rounded up
        P = get(ax,'pos');
        delete(ax);

        % Plot fluorescence w/ nuclear marker on the same plot
        [ax,p1,p2] = plotyy(X,B,X,C,'plot');
        %grid(ax(1),'on')
        %p1.LineWidth = 2;
        %p2.LineWidth = 2;
        %p2.LineSTyle = '--';
        ylabel(ax(1),'GFP');
        ylabel(ax(2),'iRFP (nuclear)');
        xlabel(ax(2),'Time in frames');
        cell_title = ['Position ',char(pos_ID),'; Cell # ',char(cell_no),'; Lifespan: ',num2str(lifespan)];
        title(cell_title);

        set(ax,'pos',P)
    end


    % Figure containing composite of traces for all cells
    %figure;










end
