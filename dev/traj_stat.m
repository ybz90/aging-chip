function traj_stat(traj1,flu1,traj2,flu2,w)

    % DESCRIPTION INPUTS ARE TWO CELL ARRAYS CONTAINING MATRICES IN FORMAT OF TRAJ_EXPORT; THEY SHOULD HAVE ONLY ONE CHANNEL EACH BUT ITS OK IF THEY HAVE MORE, ITS JUST THAT ONLY THE FIRST WILL BE USED, THOUGH I CAN ALLOW TWO MORE INPUTS TO SPECIFY THE CHANNEL TO BE USED AS WELL

    % Yuan Zhao 06/23/2015

    % TO DO: IMPLEMENT AVERAGING BY WINDOW SIZE w
    % ISSUE: WE SHOULD NOT COUNT the bud where a cell goes bye bye and buds out as the last cycle average becomes NaN for some cuz theres nothing in it.... icky icky... though i suppose it shouldnt be an issue going fwd with the chip that buds both ways. for now, made some temp lame code averaging the last two cell cycles if the last cell cycle count == X(end)?

    % TO DO 2: PUT TRAJ1, TRAJ2 INTO CELL ARRAY; USE FOR LOOP FOR THE POPULATION OF ARRAYS WHERE l=1:2, traj = trajstore{l}, and raw, fold, absolute {l,#} goes in OKAY DOOD, I GUESS WE NEED TO DO THIS FOR FLUCH TOO


    % Init cell arrays for storing values of each trajectory at start, mid, max within 3 cycles of end, and end
    %num_flu = numel(flu_array);
    raw = cell(2,4); % store raw values
    fold = cell(2,3); % store relative/fold change
    absolute = cell(2,3); % store absolute change


    % Store traj1, traj2 and flu1, flu2 in cell arrays
    all_traj = {traj1,traj2};
    all_flu = {flu1,flu2};


    % Populate raw, fold, and absolute cell arrays

    % For each set of trajectories
    for j = 1:2

        % Current set of trajectories input and flu channel to check
        curr_traj = all_traj{j};
        curr_flu = all_flu{j};

        % Total number of cells to be plotted
        num_cell_all = numel(curr_traj);

        % For every cell in curr_traj...
        for k = 1:num_cell_all

            % Get curr_cell data from curr_traj
            curr_cell = curr_traj{k};

            % Use fluoresence data for channel curr_flu in curr_traj
            flu_data = curr_cell(:,curr_flu+1);

            % Get raw values of initial, mid, max near end, and end and store in raw{l,1-4} respectively
            % Also, calculate fold and absolute change and store in fold/absolute{l,1-3} respectively

            % init
            initval = flu_data(1);
            raw{j,1} = horzcat(raw{j,1},initval);

            % get budding times
            budvals = curr_cell(:,end-1);
            X = curr_cell(:,1);
            cycles = []; % get list of budding cell times
            for m = 1:numel(budvals)
                if budvals(m) == 1
                    cycles(end+1) = X(m);
                end
            end
            cycles = sort(cycles);

            % mid
            %midpoint = ceil(length(flu_data)/2);
            %midval = flu_data(midpoint);
            % average of middle cycle
            midpoint = floor(numel(cycles)/2);
            mid_start = cycles(midpoint);
            mid_end = cycles(midpoint+1);
            midval = mean(flu_data(mid_start-X(1)+1:mid_end-X(1)+1));
            raw{j,2} = horzcat(raw{j,2},midval); % store raw value
            fold{j,1} = horzcat(fold{j,1},midval/initval); % divide by initval
            absolute{j,1} = horzcat(absolute{j,1},midval-initval); % subtract initval

            % max_end (max near end within a certain # of cell cycles)
            allow_cycles = 3; % # of cell cycles to track back from end to find max
            tr_start = cycles(end-allow_cycles+1);
            max_end = max(flu_data(tr_start-X(1)+1:end));
            raw{j,3} = horzcat(raw{j,3},max_end);
            fold{j,2} = horzcat(fold{j,2},max_end/initval);
            absolute{j,2} = horzcat(absolute{j,2},max_end-initval);

            % end
            %endval = flu_data(end);
            % average of last cell cycle
            % TEMPORARY MEASURE DUE TO RECORDING TIME OF LAST BUD == X(end) FROM BUDDING OUT
            % FOR NEW DATA WITH BI-DIRECTIONAL BUDDING, WE CAN JUST DO LAST=CYCLES(END); FLUDATA(LAST-X+1:END)
            if cycles(end) == X(end)
                last_start = cycles(end-1);
                endval = mean(flu_data(last_start-X(1)+1:end));
            else
                last_start = cycles(end);
                endval = mean(flu_data(last_start-X(1)+1:end));
            end
            raw{j,4} = horzcat(raw{j,4},endval);
            fold{j,3} = horzcat(fold{j,3},endval/initval);
            absolute{j,3} = horzcat(absolute{j,3},endval-initval);

        end

    end


    % % Plot violin plots of distributions for the three cell arrays

    % figure;
    % %title

    % % Raw values subplot
    % subplot(3,6,1);
    % % sub title
    % % legend
    % % y-axis, x-axis
    % % calculate mu, std for every matrix in cell array raw
    % mu_std_raw = cell(1,2); % holds the means and standard deviations of every mat in raw; there are up to 2 rows of matricies for each flu channel and two columns per row for mean and std respectively
    % for n = 1:4
    %     raw{1,n};
    %     mu = mean(raw{1,n});
    %     stddev = std(raw{1,n});
    %     mu_std_raw{1,1} = horzcat(mu_std_raw{1,1},mu);
    %     mu_std_raw{1,2} = horzcat(mu_std_raw{1,2},stddev);
    % end

    % % PLOT VIOLIN FOR 1 VS 2
    % a = sort(raw{1,2});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % xlim([0 8000]);

    % subplot(3,6,2);
    % a = sort(raw{2,2});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([0 8000]);

    % subplot(3,6,3);
    % a = sort(raw{1,3});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([0 8000]);

    % subplot(3,6,4);
    % a = sort(raw{2,3});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([0 8000]);

    % subplot(3,6,5);
    % a = sort(raw{1,4});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([0 8000]);

    % subplot(3,6,6);
    % a = sort(raw{2,4});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([0 8000]);


    % subplot(3,6,7);
    % a = sort(fold{1,1});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([0 14]);

    % subplot(3,6,8);
    % a = sort(fold{2,1});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([0 14]);

    % subplot(3,6,9);
    % a = sort(fold{1,2});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([0 14]);

    % subplot(3,6,10);
    % a = sort(fold{2,2});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([0 14]);

    % subplot(3,6,11);
    % a = sort(fold{1,3});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([0 14]);

    % subplot(3,6,12);
    % a = sort(fold{2,3});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([0 14]);


    % subplot(3,6,13);
    % a = sort(absolute{1,1});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([-3000 7000]);

    % subplot(3,6,14);
    % a = sort(absolute{2,1});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([-3000 7000]);

    % subplot(3,6,15);
    % a = sort(absolute{1,2});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([-3000 7000]);

    % subplot(3,6,16);
    % a = sort(absolute{2,2});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([-3000 7000]);

    % subplot(3,6,17);
    % a = sort(absolute{1,3});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','k')
    % set(gca,'view',[-90 90])
    % b = get(gca,'xlim')
    % xlim([-3000 7000]);

    % subplot(3,6,18);
    % a = sort(absolute{2,3});
    % %a = a(1:end-1);
    % histogram(a,25,'FaceColor','r')
    % set(gca,'view',[90 -90])
    % set(gca,'XTickLabel','')
    % xlim([-3000 7000]);


    % Custom colors from rgb database
    green = [0 0.5 0];
    forestgreen = [0.1328 0.5430 0.1328];
    firebrick = [0.6953 0.1328 0.1328];
    crimson = [0.8594 0.0781 0.2344];
    rgbgray = [0.5 0.5 0.5];
    dimgray = [0.4102 0.4102 0.4102];

    figure;

    subplot(3,3,1);
    a = sort(raw{1,2});
    b = sort(raw{2,2});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = 100;
    % h1 = histfit(a,25,'kernel');
    % set(h1(1),'FaceColor',green); set(h1(2),'LineStyle','--','Color','k');
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = 100;
    % h2 = histfit(b,25,'kernel');
    % set(h2(1),'FaceColor',crimson); set(h2(2),'LineStyle','--','Color','k');
    hold on
    title('Raw values (first cycle)');
    xlabel('Fluorescence (au)');
    ylabel('Frequency');
    legend('NTS1','NTS2');
    %legend('NTS2','NTS2-histfit','URA3','URA3-histfit');

    subplot(3,3,2);
    a = sort(raw{1,3});
    b = sort(raw{2,3});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = 100;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = 100;
    title('Raw values (peak near end)');
    xlabel('Fluorescence (au)');
    ylabel('Frequency');

    subplot(3,3,3);
    a = sort(raw{1,4});
    b = sort(raw{2,4});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = 100;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = 100;
    title('Raw values (last cycle)');
    xlabel('Fluorescence (au)');
    ylabel('Frequency');

    subplot(3,3,4);
    a = sort(fold{1,1});
    b = sort(fold{2,1});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = .1;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = .1;
    hold on
    y_lim = get(gca,'ylim');
    plot([1 1],y_lim,'k:');
    title('Fold change (first cycle)');
    xlabel('Fold change');
    ylabel('Frequency');

    subplot(3,3,5);
    a = sort(fold{1,2});
    b = sort(fold{2,2});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = .2;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = .2;
    hold on
    y_lim = get(gca,'ylim');
    plot([1 1],y_lim,'k:');
    title('Fold change (peak near end)');
    xlabel('Fold change');
    ylabel('Frequency');

    subplot(3,3,6);
    a = sort(fold{1,3});
    b = sort(fold{2,3});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = .2;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = .2;
    hold on
    y_lim = get(gca,'ylim');
    plot([1 1],y_lim,'k:');
    title('Fold change (last cycle)');
    xlabel('Fold change');
    ylabel('Frequency');

    subplot(3,3,7);
    a = sort(absolute{1,1});
    b = sort(absolute{2,1});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = 100;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = 100;
    title('Absolute change (first cycle)');
    xlabel('Fluorescence (au)');
    ylabel('Frequency');

    subplot(3,3,8);
    a = sort(absolute{1,2});
    b = sort(absolute{2,2});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = 100;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = 100;
    title('Absolute change (peak near end)');
    xlabel('Fluorescence (au)');
    ylabel('Frequency');

    subplot(3,3,9);
    a = sort(absolute{1,3});
    b = sort(absolute{2,3});
    %a = a(1:end-1);
    h1 = histogram(a,25,'FaceColor','k');
    h1.Normalization = 'probability';
    h1.BinWidth = 100;
    hold on
    h2 = histogram(b,25,'FaceColor','r');
    h2.Normalization = 'probability';
    h2.BinWidth = 100;
    title('Absolute change (last cycle)');
    xlabel('Fluorescence (au)');
    ylabel('Frequency');



% figure;

%     s1 = subplot(1,3,1);
%     a = sort(fold{1,1});
%     a = a(1:end-1);
%     b = sort(fold{2,1});
%     %b = b(1:end-1);
%     bihist(a,b,50);
%     % PLOT BOX PLOTS? MAYBE ADD THIS TO BIHIST
%     hold on
%     plot(xlim/2,[1 1],'k:')
%     axis off
%     p1 = get(s1,'pos');
%     p1(3) = 0.25;
%     %p1(1) = p1(1);
%     set(s1,'pos',p1);

%     s2 = subplot(1,3,2);
%     a = sort(fold{1,2});
%     a = a(1:end-1);
%     b = sort(fold{2,2});
%     %b = b(1:end-1);
%     bihist(a,b,50);
%     % PLOT BOX PLOTS? MAYBE ADD THIS TO BIHIST
%     hold on
%     plot(xlim/2,[1 1],'k:')
%     axis off
%     p2 = get(s2,'pos');
%     p2(3) = 0.25;
%     %p2(1) = p2(1) - 0.03;
%     p2(1) = p1(1) + 0.25;
%     set(s2,'pos',p2);

%     s3 = subplot(1,3,3);
%     a = sort(fold{1,3});
%     a = a(1:end-1);
%     b = sort(fold{2,3});
%     %b = b(1:end-1);
%     bihist(a,b,50);
%     % PLOT BOX PLOTS? MAYBE ADD THIS TO BIHIST
%     hold on
%     plot(xlim/2,[1 1],'k:')
%     %axis off
%     p3 = get(s3,'pos');
%     p3(3) = 0.25;
%     %p3(1) = p3(1)-0.06;
%     p3(1) = p1(1) + 0.5;
%     set(s3,'pos',p3);


%     title('NTS1 v URA3; fold change at mid avg, max end, and max avg');


    % % Fold change subplot
    % subplot(1,6,2);
    % % sub title
    % % legend
    % % y-axis, x-axis
    % % calculate mu, std for every matrix in cell array fold

    % % PLOT VIOLIN FOR 1 OR 2

    % % Raw values subplot
    % subplot(1,6,3);
    % % sub title
    % % legend
    % % y-axis, x-axis
    % % calculate mu, std for every matrix in cell array absolute


    % Plot scatterplot correlating lifespan with values, fold change, and absolute change
    %%%%

end
