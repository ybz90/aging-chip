function export_xlsx(pos_str,all_traj,all_lifespan,flu_array)

    % export_xls(pos,all_traj,all_lifespan,[1 2 3])
    % duplicate lifespan.txt, remove the ones w/ ERC and recomb or are just really low bg
    % export FULL traces for the rest

    % Yuan Zhao 05/18/2015

    fluN = 3;

    % one giant cell array with horz concatenated trajectories for all the cells frames is rows, cols is cell, and each cell mat is a flu ch)
    traj_export = cell(1,fluN);
    for p = 1:fluN
        ars = 0:900;
        %ars = zeros(1,901);
        ars = ars';
        traj_export{p} = ars;
    end

    %array for storing cycles
    max_cycles = cell(1,numel(all_lifespan)); % the highest # of rep life cycles for each position to be plotted; simply take the # of bud frames as the traps with fewer than the max will fill in a 0 anyway
    min_cycles = cell(1,numel(all_lifespan));
    for c = 1:numel(pos_str)
        curr_life = all_lifespan{c};
        max_cycles{c} = numel(curr_life(1,5:end));
    end
    max_cycles = max(cell2mat(max_cycles))

    %life_export = zeros(1,max_cycles+1);
    life_export = 0:max_cycles;


    count = 1;

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

            for z = 1:numel(flu_array)

                a = curr_traj(1:900,cell_ID,flu_array(z));
                %b = [str2num(pos); cell_ID];
                b = count;
                a = vertcat(b,a);
                %pos, cell_ID, size(a)

                c = flu_array(z);

                traj_export{c} = horzcat(traj_export{c},a);

                % SAVE IN RESPECTIVE HORZCAT TO ARRAYS BY FLU CH, add to front the pos and cell_ID
            end

            d = curr_life(l,5:end);
            missing = zeros(1,max_cycles-numel(d));
            d = horzcat(count,d,missing);
            life_export = vertcat(life_export,d);



            count = count+1;

        end

    end

    %size(traj_export)
    %traj_export
    %traj_export{1}
    %traj_export{2}
    %traj_export{3}

    %life_export

    xlswrite('20150521_NTS2-NTS2-c2_mCherry.xlsx',traj_export{1}','C2 - mCherry');
    xlswrite('20150521_NTS2-NTS2-c3_GFP.xlsx',traj_export{2}','C3 - GFP');
    %xlswrite('20150521_NTS1-NTS2-c4.xlsx',traj_export{3},'C4 - Cy5');

    xlswrite('20150521_NTS2-NTS2-bud_times.xlsx',life_export,'Budding Times');
end
