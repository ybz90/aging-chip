function traj_stat(pos_str,all_traj,all_lifespan,flu_array)

    % DESCRIPTION

    % Yuan Zhao 06/03/2015

    fluN = 3;

    %array for storing cycles
    num_traj = 0;
    chrono_len = [] %chronological lifespans
    for c = 1:numel(pos_str)
        curr_life = all_lifespan{c}
        for d = 1:numel(curr_life(1:end,1))
            chrono_len = horzcat(chrono_len,[curr_life(d,4)-curr_life(d,3)]);
            num_traj = num_traj + 1;
        end
    end
    chrono_len
    shortest_chrono = min(chrono_len)
    num_traj


    % all shortened paths from death backwards
    short_traj = cell(1,num_traj);

    count = 1
    for k = 1:numel(pos_str)
        pos = pos_str{k}; % set current position
        curr_traj = all_traj{k}; %load current trajectory file based on position index
        curr_life = all_lifespan{k}; %load lifespan data

        sz = size(curr_life);
        num_cells = sz(1); % number of rows in current position's lifespan data, ie # of cell trajectories to plot
        for l = 1:num_cells
            cell_ID = curr_life(l,1);
            curr_dead = curr_life(l,4)-shortest_chrono:curr_life(l,4);
            a = curr_traj(curr_dead,cell_ID,1);
            short_traj{count} = a;
            count = count + 1;
        end
    end


    %short_traj{1}

    % first teh quick and dirty of the ugly stats ie one color analysis

    %define window size tau
    %correct by subtraction of of window avg of end point; move all curves up or down to this avg
    %window avg of substraced starts
    %find std dev of ensemble avg of window start avgs
    %check start avg + 2 overall ensemble std vs final avg; set this actually as a thrsh param

    % do this for 1, 2, both? compare efficiency



end
