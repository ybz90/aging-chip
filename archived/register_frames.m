function register_frames(pos,imN)

    % register_frames.m is used to register the positions of each trap across the frames to account for movement in the stage in the form of shaking or drifting over time.
    % This function is an OPTIONAL step before mask generation and cell trajectory tracking and need not be used if the stage is stable throughout the course of the experiment.

    % Original code by Meng Jin; last edited by Yuan Zhao 05/05/15


    for imid = 1:imN

        fprintf('Registering frame number %d.\n', imid); %debug

        % Load the names and output dirs of the raw frames to be registered, for the current time frame (imid)
        % write output to the same path as input, overwriting raw images
        thrsh_name = ['xy',pos,'/c1_thr/xy',pos,'_c1_thr_t',sprintf('%04g',imid),'.tif'];
        I0 = imread(thrsh_name);

        phase_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        Iph = imread(phase_name);

        fluor_name = ['xy',pos,'/c2/xy',pos,'_c2_t',sprintf('%04g',imid),'.tif'];
        Iflr = imread(fluor_name);

        nuclr_name = ['xy',pos,'/c3/xy',pos,'_c3_t',sprintf('%04g',imid),'.tif'];
        Inuc = imread(nuclr_name);


        I1 = imopen(I0, strel('rectangle', [2, 50])); %changed from [2, 100]; try [2, 75]?
        I1_extend = imdilate(I1, strel('line',500,0));

        % Label connected components in dilated image; after removing small objects
        I1_lb = bwlabel(bwareaopen((I1_extend==0),5000));
        center_blck = (I1_lb==2);

        I1_cnt = (I1_lb==2).*(I0>0);
        I2 = imopen(I1_cnt, strel('rectangle',[10, 2]));
        I3 = imdilate(I2, strel('rectangle',[100, 15]));
        Ia = I3.*I1_cnt;

        Ia1 = imdilate(Ia, strel('disk',6));
        Ia2 = imfill(Ia1, 'holes');
        % Remove any non-column features
        Ia3 = bwareaopen(Ia2, 1000);

        [Ia_lb,colN]= bwlabeln(Ia3);

        center_xy = regionprops(center_blck, 'Centroid');
        column_xy = regionprops(Ia_lb, 'Centroid');

        columnN = length(column_xy);
        curr_xy = round([column_xy(1).Centroid(1), center_xy.Centroid(2)]);

        % Initialize first frame, set initial position
        if imid == 1
            use_refxy = curr_xy;

            Imask_new = I0;
            Iflr_new = Iflr;
            Inuc_new = Inuc;
            Iph_new  = Iph;

        % Update positions and shift flu, nuc, phase, and threshold images to align
        else
            d_xy = curr_xy - use_refxy;
            d_col = d_xy(1);
            d_row = d_xy(2);

            Iph_new = uint16(zeros(size(I0)));
            Iflr_new = Iph_new;
            Inuc_new = Iph_new;

            Imask_new = uint8(zeros(size(I0)));

            COL = size(I0,2);
            ROW = size(I0,1);

            if d_col >0
                new_col= 1:COL - d_col;
                curr_col= d_col+1:COL;
            else
                new_col= -d_col+1:COL;
                curr_col= 1:COL + d_col;
            end

            if d_row >0
                 new_row= 1:ROW - d_row;
                 curr_row = d_row+1:ROW;
            else
                new_row = -d_row+1:ROW;
                curr_row =1: ROW + d_row;
            end

            Imask_new(new_row, new_col) = I0(curr_row, curr_col);
            Iflr_new(new_row, new_col) = Iflr(curr_row, curr_col);
            Inuc_new(new_row, new_col) = Inuc(curr_row, curr_col);
            Iph_new(new_row, new_col) = Iph(curr_row, curr_col);
        end

        % Export new, positionally registered images
        imwrite(Imask_new, thrsh_name);
        imwrite(Iflr_new, fluor_name);
        imwrite(Inuc_new, nuclr_name);
        imwrite(Iph_new, phase_name);
    end
end
