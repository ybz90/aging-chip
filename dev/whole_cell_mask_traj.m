function whole_cell_mask_traj(pos,imN,fluN)

    % mask_traj.m is used to generate a mask from the nuclear marker images and produce trajectories for the mother cells in each trap.

    % Yuan Zhao 05/19/2015



    % http://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
    % INITIALIZATION %
    colN = 7;

    % Create output directory for storing mother mask+phase overlays
    mkdir(strcat('xy',pos,'/mother_overlay'))

    % Get image dimensions
    ph_1 = ['xy',pos,'/c1/xy',pos,'_c1_t0001.tif'];
    I_ph_1 =imread(ph_1);
    sz = size(I_ph_1); % rows (image height), cols (image width)
    height = sz(1);
    width = sz(2);

    % Define column slice width as the image width divided by colN
    block = round(width/colN);

    % Store the centroids and column masks of each mother cell
    for w = 1:colN
        mother_x(w) = block/2; %initialize with x-centroid right in the middle of block
    end
    mother_y = zeros(1,colN);
    mother_BW = cell(1,colN); %initalize with completely blank columns
    for p = 1:colN
        mother_BW{p} = zeros(height,block);
    end


    % FOR EACH FRAME... %
    for imid = 1:imN

    % STEP 1: MASK GENERATION %
        fprintf('Generating masks and trajectories for xy%s, frame %d.\n', pos, imid); %debug

        % Input directory and image paths for nuclear marker, fluorescence, and phase channel images
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph_0 = imread(ph_name);

        % The nuclear marker channel is, by convention, the last fluorescent channel
        nuc_name = ['xy',pos,'/c',num2str(fluN+1),'/xy',pos,'_c',num2str(fluN+1),'_t',sprintf('%04g',imid),'.tif'];
        I_nuc_0 = imread(nuc_name);

        % flu c2
        flu_name = ['xy',pos,'/c2/xy',pos,'_c2_t',sprintf('%04g',imid),'.tif'];
        %flu_name = ['xy',pos,'/c3/xy',pos,'_c3_t',sprintf('%04g',imid),'.tif'];
        I_flu_0 = imread(flu_name);

        % Output directory and mask image name
        out_name_mother = ['xy',pos,'/mother_overlay/xy',pos,'_mother_overlay_t',sprintf('%04g',imid),'.tif'];

        % Initialize mask image for this frame
        if rem(width,colN) == 0 % if image width is divisible by 7, init blank mask array
            I_nuc_mask = [];
            I_mother_mask = [];
        else
            I_nuc_mask = zeros(height,1); % since 512 is not divisible by 7, the block size (col width) is rounded and we need to add a column of black pixels to get 511+1 = 512 px
            I_mother_mask = zeros(height,1,3); %
        end


        % FOR EACH COLUMN IN THE CURRENT FRAME... %
        for i = 1:colN
            I_nuc = I_nuc_0(:,1+(i-1)*block:i*block); %current column of nuclear marker image
            %figure; imagesc(I_nuc(:,1+(i-1)*block:i*block)) %debug
            I_ph = I_ph_0(:,1+(i-1)*block:i*block);


% PHASE CELL OUTLINES + FEATURES FILLED IN
            % make a phase mask to subtract from flu
            %[levelph EMph] = graythresh(I_ph);
            bwp = im2bw(I_ph,0.0275+.011);
            %bwp = imclose(bwp,strel('disk',2));
            %figure; imshow(bwp);

            % take initial phase image... and fill in the trap area completely; erode it to roughly the width of the actual trap, then invert it so the features are filled in and not the trap... then add it to the initial phase to basically fill in the features, leaving only the cells and a few artifacts not filled in
            % this will act as a boundary, a cookie mold, for the fluorescence mask to constrain it
            % then we will watershed the phase-constrained fluorescence image, with n uclear markers as marker tags

            %trap = imdilate(bwp,strel('line',10,90));
            trap = imdilate(bwp,strel('rectangle',[10,15]));
            trap = imdilate(trap,strel('disk',5));
            trap = imfill(trap,'holes');
            trap = imerode(trap,strel('disk',12));
            trap = imcomplement(trap);
            %figure; imshow(trap)

            bwp = bwp + trap;
            bwp = im2bw(bwp);
            %figure; imshow(bwp)
            bwp = imclose(bwp,strel('disk',2));
            %figure; imshow(bwp)


% FLU
            I_flu = I_flu_0(:,1+(i-1)*block:i*block);
            %figure; imshow(I_flu)

            % background = imopen(I_flu,strel('disk',50));
            % figure
            % surf(double(background(1:8:end,1:8:end))),zlim([0 255]);
            % set(gca,'ydir','reverse');


            % I_flu = imadjust(I_flu);
            % %figure; imshow(I_flu);

            % % I_FLU HIST
            % [maxval maxloc] = max(I_flu(:));
            % [maxloc_row maxloc_col] = ind2sub(size(I_flu), maxloc);
            % maxval;

            % %medval = median(I_flu(:))
            % I_flu_list = I_flu(:);
            % I_flu_list = sort(I_flu_list);
            % ind_num = numel(I_flu_list); %number of pixels in each column
            % I_flu_list(floor(ind_num*0.5):end);

            % %y = I_flu_list(floor(ind_num*0.99))
            % z = I_flu_list(floor(ind_num*0.8))
            % % 1533 - 113; level of 0.0333... so, lets try (max-0.8point)/(1533-113) * 0.0333?


              I_flu = adapthisteq(I_flu);
              [level EM] = graythresh(I_flu);


              % % adfasdf
              % ref = 1533-113;
              % level2 = (ref-double(maxval-z))/ref*0.0333
              % level = level+level2


              %level = 0.2;
              bw = im2bw(I_flu,level);
              bw2 = imfill(bw,'holes');
              bw3 = imdilate(bw2, strel('disk',1)); %dilate mask with disk
              bw4 = bwareaopen(bw3, 50);
              bw4_perim = bwperim(bw4);
              overlay1 = imfuse(I_flu,bw4_perim,'falsecolor','ColorChannels','red-cyan');

              %figure; imshow(BW3) %debug
              %figure; imshow(overlay1)
             %% figure; imshow(bw4);

              %bw4 = imerode(bw4,strel('disk',1));
              %bw4 = bwareaopen(bw4,25);
              %bw4 = imerode(bw4,strel('disk',2));
              %bw4 = bwareaopen(bw4,35);
              %figure; imshow(bw4)

% SUBTRACT PHASE FROM FLU TO GET ONLY THE AREAS WITH FLU THRESH WITHIN THE BOUNDS OF PHASE
              bwsub = bw4-bwp;
              bwsub = bwsub>0;
             %bwsub = im2bw(bwsub);
              %figure; imshow(bwsub);

% FIND NUC MARKER
            Icf = I_nuc; %current column of nuclear marker image
            %figure; imagesc(I_nuc(:,1+(i-1)*block:i*block)) %debug

            % INITIAL THRESHOLDING %
            % Otsu's threshold nuclear marker image to obtain binary column mask
            %[level EM] = graythresh(Icf);
            level = 0.05;
            BW = im2bw(Icf,level);
            BW2 = imfill(BW,'holes');
            BW3 = imdilate(BW2, strel('disk',1)); %dilate mask with disk
            %figure; imshow(BW3) %debug

            mask_prop = regionprops(BW3,Icf,'Area','Centroid'); %areas and centroids of the cell masks for the current column

            % Get centroids of all cells and their x,y coordinates
            all_centroids = [mask_prop(:).Centroid];
            x_centroids = all_centroids(1:2:end-1);
            y_centroids = all_centroids(2:2:end);

            % Structure for holding the cell fluorescence and other property data for all cells
            % mask_prop_2 has dimensions of 3 rows (cell properties) x n cols (# of cells in the trap)
            mask_prop_2 = [mask_prop(:).Area; x_centroids; y_centroids];
            %fprintf('Column #%d has %d cells.\n', i, length(mask_prop_2(:))/4); %debug

            % The centroid with the highest y-coordinate value is the cell that is lowest in the trap, aka the mother cell; its index, idx, specifies the column in mask_prop_2 that stores its properties
            [max_y,idx] = max(y_centroids);

            % Properties of the mother cell, from the column #idx
            mother_prop = mask_prop_2(:,idx); %mother_prop(1) is the area of the mother cell
            areas = mask_prop_2(1,:);


            % MASK OPTIMIZATION TO IDENTIFY HARD-TO-DETECT MOTHER CELLS %
            % If there are no cells detected via nuclear mask in this column, reduce threshold up to 3 times to try to detect them
            num_tries = 0;
            while isempty(mother_prop) && num_tries < 3
                num_tries = num_tries + 1;
                level = level-0.005;
                BW = im2bw(Icf,level); %repeat Otsu's method with new paramaters
                BW2 = imfill(BW,'holes');
                BW3 = imdilate(BW2, strel('disk',1));
                mask_prop = regionprops(BW3,Icf,'Area','Centroid');

                all_centroids = [mask_prop(:).Centroid];
                x_centroids = all_centroids(1:2:end-1);
                y_centroids = all_centroids(2:2:end);
                mask_prop_2 = [mask_prop(:).Area; x_centroids; y_centroids];
                [max_y,idx] = max(y_centroids);
                mother_prop = mask_prop_2(:,idx);
                areas = mask_prop_2(1,:);
            end

            % Check the following conditions to determine whether or not to further reduce the thresold level, thereby increasing the mask radius
            % NOTE: Essentially, we want to maximize the likelihood of capturing a mother cell mask while avoiding oversaturation, resulting in background artifacts. The mother cell only needs to be larger than the areaopen threshold in the next step.
            while ~isempty(mother_prop) && (mother_prop(1) < 35 | min(areas) < 25) && max(areas) < 150 && num_tries < 8
                % Check that there is a mother cell identified
                % Check if the mother cell is too small, OR if any other cell is also too small (this or statement deals with the scenario wherein the lowermost cell of the initial threshold may be sufficiently large to skip the while loop, but it is not the actual mother cell, which may not be detected that that thrsh level)
                % Stop if the largest cell exceeds too large a size, as this could imply oversaturation during threshold
                % Limit the number of retries to < 3
                level = level-0.005; %reduce level for Otsu's method
                num_tries = num_tries + 1;
                BW = im2bw(Icf,level); %repeat Otsu's method with new paramaters
                BW2 = imfill(BW,'holes');
                BW3 = imdilate(BW2, strel('disk',1));
                mask_prop = regionprops(BW3,Icf,'Area','Centroid');

                all_centroids = [mask_prop(:).Centroid];
                x_centroids = all_centroids(1:2:end-1);
                y_centroids = all_centroids(2:2:end);
                mask_prop_2 = [mask_prop(:).Area; x_centroids; y_centroids];
                [max_y,idx] = max(y_centroids);
                mother_prop = mask_prop_2(:,idx);
                areas = mask_prop_2(1,:);

                %[mother_prop(1),min(areas),max(areas),num_tries] %debug
            end
              % MOTHER
              %%figure; imshow(BW3)

              % % AUTO GEN W/O NUC MARKER?
              % BW3 = imextendedmin(D,2);



% ADD NUCLEAR MARKER TO THE BWSUB MASK SO FAR TOO IN CASE ITS NOT INCLUDED
              bwsub = bwsub + BW3;
              bwsub = im2bw(bwsub);
              %figure; imshow(bwsub)
              bwsub = double(bwsub);
              %figure; imshow(bwsub)


              A = imimposemin(bwsub,BW3);
              %figure; imshow(A); %test show min

% THRESHOLD THE SUBTRACTED FLU IAMGE WITH MOTHER CELL AS MARKER
              %D = -bwdist(~bw4);
              D = -bwdist(~bwsub);
              %figure; imshow(D,[])
              % Ld = watershed(D);
              % figure; imshow(label2rgb(Ld));

              %%% figure; imshowpair(bw4,BW3,'blend')

              %figure; imshowpair(bwsub,BW3,'blend')

% V1
              D2 = imimposemin(D,BW3);
              Ld2 = watershed(D2);
              bwZ = bwsub;
              bwZ(Ld2 == 0) = 0;
              bwZ = bwareaopen(bwZ,75);
              %figure; imshow(bwZ)
% END V1

              % % TAKE MOTHER ONLY AFTER WATERSHED SEGMENT
              % bwZ = bwZ';
              % [P,num] = bwlabel(bwZ);
              % bwZ = (P==num);
              % bwZ = bwZ';

              % bwZ = imdilate(bwZ,strel('disk',2));


              %overlay with phase
              bwZperim = bwperim(bwZ);
              final = imfuse(I_ph,bwZperim,'falsecolor','ColorChannels','red-cyan');
              %figure; imshow(final)

              %%size(final)
              %%sz
              % Add current mother cell mask to the overall mother mask image
              I_mother_mask = horzcat(I_mother_mask,final);






            % % MINIMUM SPANNING TREES METHOD
            % pred = 4; % use nuc markers to find # of marker objects
            % [Tree, pred] = graphminspantree(I_ph)



            % [~, threshold] = edge(I, 'sobel');
            % fudgeFactor = .5;
            % BWs = edge(I,'sobel', threshold * fudgeFactor);
            % figure, imshow(BWs), title('binary gradient mask');

            % se90 = strel('line', 3, 90);
            % se0 = strel('line', 3, 0);

            % BWsdil = imdilate(BWs, [se90 se0]);
            % figure, imshow(BWsdil), title('dilated gradient mask');

            % BWdfill = imfill(BWsdil, 'holes');
            % figure, imshow(BWdfill);
            % title('binary image with filled holes');

            % BWnobord = imclearborder(BWdfill, 4);
            % figure, imshow(BWnobord), title('cleared border image');

            % seD = strel('diamond',1);
            % BWfinal = imerode(BWnobord,seD);
            % BWfinal = imerode(BWfinal,seD);
            % figure, imshow(BWfinal), title('segmented image');


            end

        % Output mother mask+phase overlay
        imwrite(I_mother_mask, out_name_mother)
        %%figure; imshow(I_mother_mask)
    end

end
