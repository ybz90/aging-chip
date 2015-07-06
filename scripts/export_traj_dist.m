function export_traj_dist(curr_traj,curr_flu,w)

    % INPUT IS A CELL ARRAY OF STMs (IDEALLY AFTER COMBINING ALL RELEVANT DATA TOGETHER WITH combine_traj S.T. EACH STM HAS ONLY ONE FLU CHANNEL e.g. NTS2 all data; ie traj_all is input).
    % OUTPUT is distributions of first, 1/4, mid, 3/4, and end cycles (per w window size); export as csv with a csv for each distribution.
    % FOR NOW, MANUALLY SPECIFY WHETHER TO OUTPUT RAW, FOLD CHANGE, OR ABS CHANGE DATA.

    % TO DO: IMPLEMENT AVERAGING BY WINDOW SIZE w
    % ISSUE: WE SHOULD NOT COUNT the bud where a cell goes bye bye and buds out as the last cycle average becomes NaN for some cuz theres nothing in it.... icky icky... though i suppose it shouldnt be an issue going fwd with the chip that buds both ways. for now, made some temp lame code averaging the last two cell cycles if the last cell cycle count == X(end)?


    % Init cell arrays for storing values of each trajectory at (1st cycle, .25*nth cycle, mid cycle, 0.75*nth cycle, nth cycle)
    %num_flu = numel(flu_array);
    raw = cell(1,5); % store raw values
    fold = cell(1,5); % store relative/fold change
    absolute = cell(1,5); % store absolute change


    % Populate raw, fold, and absolute cell arrays

    % Total number of cells to be plotted
    num_cell_all = numel(curr_traj);

    % For every cell in curr_traj...
    for k = 1:num_cell_all

        % Get curr_cell data from curr_traj
        curr_cell = curr_traj{k};

        % Use fluoresence data for channel curr_flu in curr_traj
        flu_data = curr_cell(:,curr_flu+1);

        % get indices of budding times
        budvals = curr_cell(:,end-1);
        X = curr_cell(:,1); %time frames
        cycles = []; % get list of budding cell time indices
        for m = 1:numel(budvals)
            if budvals(m) == 1
                cycles(end+1) = m;
            end
        end
        cycles = sort(cycles);

        % Get raw values, fold change, and absolute change of initial, 1/4, mid, 3/4, and end and store in raw{l,1-4} and fold/absolute{l,1-3} respectively
        % init
        init_val = mean( flu_data(1:cycles(1)) );
        raw{1} = horzcat(raw{1},init_val);

        fold{1} = horzcat(fold{1},1);
        change{1} = horzcat(fold{1},1);

        % 1/4th cycle
        quart_point = floor(numel(cycles)/4);
        quart_start = cycles(quart_point);
        quart_end = cycles(quart_point+1);
        quart_val = mean( flu_data(quart_start:quart_end) );
        raw{2} = horzcat(raw{2},quart_val);
        fold{2} = horzcat(fold{2},quart_val/init_val);
        absolute{2} = horzcat(absolute{2},quart_val - init_val);

        % mid cycle
        mid_point = floor(numel(cycles)/2);
        mid_start = cycles(mid_point);
        mid_end = cycles(mid_point+1);
        mid_val = mean( flu_data(mid_start:mid_end) );
        raw{3} = horzcat(raw{3},mid_val);
        fold{3} = horzcat(fold{3},mid_val/init_val);
        absolute{3} = horzcat(absolute{3},mid_val - init_val);

        % 3/4th cycle
        thfo_point = floor(3*numel(cycles)/4);
        thfo_start = cycles(thfo_point);
        thfo_end = cycles(thfo_point+1);
        thfo_val = mean( flu_data(thfo_start:thfo_end) );
        raw{4} = horzcat(raw{4},thfo_val);
        fold{4} = horzcat(fold{4},thfo_val/init_val);
        absolute{4} = horzcat(absolute{4},thfo_val - init_val);

        % average of last cell cycle
        % TEMPORARY MEASURE DUE TO RECORDING TIME OF LAST BUD == X(end) FROM BUDDING OUT
        % FOR NEW DATA WITH BI-DIRECTIONAL BUDDING, WE CAN JUST DO LAST=CYCLES(END); FLUDATA(LAST-X+1:END)
        if cycles(end) == X(end)
            last_start = cycles(end-1);
            endval = mean(flu_data(last_start:end));
        else
            last_start = cycles(end);
            endval = mean(flu_data(last_start:end));
        end
        raw{5} = horzcat(raw{5},endval);
        fold{5} = horzcat(fold{5},endval/init_val);
        absolute{5} = horzcat(absolute{5},endval-init_val);

    end

    % export distribution of raw, fold, or abs as a separate .csv for each of the five time points
    % NOTE: for now, manually edit this code to choose which one to export and the export name

    export_name = 'raw_stdnorm';

    for r = 1:numel(raw)
        xlswrite([export_name,'_',num2str(r)],raw{r});
    end

end
