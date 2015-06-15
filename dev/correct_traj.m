function correct_traj(pos_str,all_traj,all_lifespan,flu_array)

    % function to correct concentration/intensity trajectories for dilution
    % correct_traj(pos,all_traj,all_lifespan,[1 2])


    % CONTIANERS FOR SAVING CORRECT TRAJ FOR EACH CHANNEL GO HERE
    %


    % For every position in pos_str
    for k = 1:numel(pos_str)

        pos = pos_str{k}; % set current position
        curr_traj = all_traj{k}; %load current trajectory file based on position index
        curr_life = all_lifespan{k}; %load lifespan data

        sz = size(curr_life);
        num_cells = sz(1); % number of rows in current position's lifespan data, ie # of cell trajectories to plot


        % For every cell in the current traj to be plotted
        for l = 1:num_cells

            cell_ID = curr_life(l,1); %the trap in which the current cell is actually located, as not every cell is to be plotted and has info in curr_life


            for z = 1:numel(flu_array)

                cell_traj = curr_traj(:,cell_ID,flu_array(z));

                % SMOOTH?
                cell_traj2 = [];
                for q = 1:numel(cell_traj)-1
                    cell_traj2(end+1) = mean([cell_traj(q),cell_traj(q+1)]);
                end
                cell_traj2(end+1) = cell_traj2(end);
                cell_traj = cell_traj2;

                numel(cell_traj)


                lifespan = curr_life(l,5:end); %lifespan budding times
                % sort lifespan in ascending order and remove zeros
                lifespan = lifespan(lifespan>0);
                lifespan = sort(lifespan);
                %lifespan = lifespan-8; % DEBUG SHIFT OVER RIGHT TO CAPTURE THE FULL VALLEY?

                %lifespan(1), lifespan(2)-1
                %cell_traj(lifespan(1):lifespan(2)-1)

                % OK SO NOW WE HAVE THE TRACJ AND THE LIFESPANS, LETS DO THE CORRECTION ALGORITHM FOR EACH CHANNEL WE CARE ABOUT

                corr_cell_traj = []; % init corrected cell traj
                partition_adding = 0; % amount to add each time



                for m = 1:numel(lifespan)-1 % for each subpartition
                    part_start = lifespan(m); %
                    part_end = lifespan(m+1)-1; %

                    % current partition
                    curr_part = cell_traj(part_start:part_end);
                    %curr_part2 = cell_traj2(part_start:part_end);


                    %find min/max and their indices (aka time)
                    [part_max max_ind] = max(curr_part);
                    %[part_min min_ind] = min(curr_part(max_ind+1:end)); % find min of the partition AFTER the max
                    % for now arbitrarily defining min_ind as the end of the partition
                    partsize = numel(curr_part);

                    if max_ind < .25*partsize | max_ind > .75*partsize
                        max_ind = floor(partsize/2);
                        part_max = curr_part(max_ind);
                    end

                    min_ind = partsize;

                    %k,l,m
                    %time from max to min
                    t_delta = min_ind - max_ind;
                    a = 1/t_delta;

                    % CELL CURR_PART UP TO AREA TO BE CORRECTED (ie DILUTION ZONE)
                    % CORRECT THIS PARTITION BTW MAX AND MIN; ADD TO CORR_CELL_TRAJ
                    for s = 1:max_ind
                        corr_cell_traj(end+1) = curr_part(s) + partition_adding;
                    end
                    for t = max_ind+1:partsize%min_ind %1:t_delta
                        t_corr = t - max_ind;
                        correction = double(-part_max* 2^(-a*t_corr) + part_max);
                        corr_cell_traj(end+1) = curr_part(t) + correction + partition_adding;
                    end
                    % % TAKE CURR_PART AFTER CORRECTED AREA + C_DELTA
                    % for u = min_ind+1:numel(curr_part)
                    %     corr_cell_traj(end+1) = curr_part(u) + partition_adding + 0.5*part_max;
                    % end

                    % % INCREASE CELL_TRAJ ACCORDINGLY
                    % cell_traj = cell_traj + c_delta;
                    partition_adding = partition_adding + 0.5*part_max;


                    % corr_cell_traj = horzcat(corr_cell_traj,curr_part);


                    %cell_traj(100:end)
                    %c_delta(100:end)



                end

                % POST CORRECTION ALGORITHM; NOW LETS SAVE THE CORRECT TRAJ TO A LARGE CELL FOR CORRECTED TRAJ PER FLU CHANNEL
                %

                % DEBUG: PLOT CORRECTED CELL TRAJ VS CELL TRAJ
                numel(corr_cell_traj)
                numel(cell_traj(lifespan(1):lifespan(end)-1))

                j = lifespan(1):lifespan(end)-1;


                % THE CONTINUOUS STEADY INCREASE CURVE; TAKE SLOPE OF FIRST CELL CYCLE BEFORE DILUTION AND JUST PLOT THAT ALL THE WAY; THIS IS WHAT IT WOULD LOOK LIKE IF THERE WAS NO UNSILENCING
                poo = cell_traj(lifespan(1):lifespan(2)-1);
                [poo_max poo_ind] = max(poo);
                slope = (max(poo)-poo(1))/(poo_ind-1);
                poo_plot = slope*(j-(lifespan(1)-1))+poo(1);

                %corr_cell_traj == cell_traj(j)
                %cell_traj(j);

                figure;
                %plotyy(j,corr_cell_traj,j,cell_traj(j));
                %plotyy(j,corr_cell_traj,j,poo_plot);
                plot(j,corr_cell_traj,'b');
                hold on
                plot(j,cell_traj(j),'g');
                hold on
                y1 = get(gca,'ylim');
                for v = 1:numel(lifespan)
                    w = lifespan(v);
                    line([w w],y1,'Color','k','LineStyle','--')
                end
                plot(j,poo_plot,'r')
            end
        end
    end


end
