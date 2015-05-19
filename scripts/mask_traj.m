function mask_traj(pos,imN,colN,fluN)

    % mask_traj.m is used to generate a mask from the nuclear marker images and produce trajectories for the mother cells in each trap.

    % Yuan Zhao 05/15/2015


    % % Create output directory for storing masks
    % mkdir(strcat('xy',pos,'/mask'))
    % % Create output directory for storing mask+phase overlays
    % mkdir(strcat('xy',pos,'/mask_overlay'))
    % Create output directory for storing mother mask+phase overlays
    mkdir(strcat('xy',pos,'/mother_overlay'))

    % Initialize empty 3D matrix to store trajectories (frame #, trap #, flu channel #)
    traj = zeros(imN, colN, fluN);
    % Output trajectory data for all traps on this xy position
    output_data =['xy',pos,'/xy',pos,'_traj.mat']; %include position and date in the filename

    % Get image dimensions
    ph_1 = ['xy',pos,'/c1/xy',pos,'_c1_t0001.tif'];
    I_ph_1 =imread(ph_1);
    I_size = size(I_ph_1); % rows (image height), cols (image width)
    height = I_size(1);
    width = I_size(2);

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


    % For each frame...
    for imid = 1:imN

        fprintf('Frame number %d has %d columns.\n', imid, colN); %debug

        % Input directory and image paths for nuclear marker, fluorescence, and phase channel images
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph = imread(ph_name);

        % The nuclear marker channel is, by convention, the last fluorescent channel
        nuc_name = ['xy',pos,'/c',num2str(fluN+1),'/xy',pos,'_c',num2str(fluN+1),'_t',sprintf('%04g',imid),'.tif'];
        I_nuc = imread(nuc_name);

        % Input directory and image paths for every fluorescent channel
        flu_names = cell(1,fluN);
        I_flu = cell(1,fluN);
        for z = 1:fluN
            flu_names{z} = ['xy',pos,'/c',num2str(1+z),'/xy',pos,'_c',num2str(1+z),'_t',sprintf('%04g',imid),'.tif'];
            I_flu{z} = imread(flu_names{z});
        end

        % Output directory and mask image name
        % out_name = ['xy',pos,'/mask/xy',pos,'_mask_t',sprintf('%04g',imid),'.tif'];
        % out_name_overlay = ['xy',pos,'/mask_overlay/xy',pos,'_mask_overlay_t',sprintf('%04g',imid),'.tif'];
        out_name_mother = ['xy',pos,'/mother_overlay/xy',pos,'_mother_overlay_t',sprintf('%04g',imid),'.tif'];

        % Initialize mask image for this frame
        if rem(width,colN) == 0 % if image width is divisible by 7, init blank mask array
            I_nuc_mask = [];
            I_mother_mask = [];
        else
            I_nuc_mask = zeros(height,1); % since 512 is not divisible by 7, the block size (col width) is rounded and we need to add a column of black pixels to get 511+1 = 512 px
            I_mother_mask = zeros(height,1); %
        end


        % For each column in the current frame...
        for i = 1:colN
            Icf = I_nuc(:,1+(i-1)*block:i*block); %current column of nuclear marker image
            %figure; imagesc(I_nuc(:,1+(i-1)*block:i*block)) %debug


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

            %num_tries = 0;
            % Check the following conditions to determine whether or not to further reduce the thresold level, thereby increasing the mask radius
            while ~isempty(mother_prop) && (mother_prop(1) < 45 | min(areas) < 35) && max(areas) < 150 && num_tries < 8
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


            % Get the mask of the mother cell only for this column
            % bwlabel does column-wise search by default; so to do row-wise searching for the lowest object, we transpose the BW3 binary image input and then transpose back the output of bwlabel
            BW4 = BW3.';
            BW4 = bwareaopen(BW4,10); % Remove tiny objects/artifacts, BUT be careful not to accidentally remove mother masks that are very small!

            %NOTE: Add code for declumping cells; since we only care about the lowest mother cell, this will only be needed for cases where a mother and a daughter are still attached when the frame was taken, and so we can approach this using a noncircularity method to identify these scenarios. Then we can do binary watershed or something similar and take the lower most object as the mother.

            [L,num] = bwlabel(BW4); %

            if num == 0 %if there are no cells in the column, use a pure black, empty column (all zeros)
                temp_mother = zeros(height,block);
            else
                temp_mother = (L==num); %image is only the mask of the mother cell
                temp_mother = temp_mother.'; %transpose to orient axes correctly
            end

            % Increase the size of the mother if it is too small
            mother_area = regionprops(temp_mother,'Area');
            if ~isempty(mother_area)
                mother_area_2 = [mother_area(1).Area];
                while mother_area_2 < 50
                    temp_mother = imdilate(temp_mother, strel('disk',1));
                    mother_area = regionprops(temp_mother,'Area');
                    mother_area_2 = [mother_area(1).Area];
                end
            end

            % Find centroid of current col mother cell, update mother_x, mother_y, and mother_BW arrays IF the current mother cell meets the following criteria:
            mother_prop = regionprops(temp_mother,'centroid'); % find centroid of the current column mother cell
            if ~isempty(mother_prop) %if there is a mother cell in this particular column
                mother_prop2 = [mother_prop(1).Centroid(1),mother_prop(1).Centroid(2)];
                % figure; imshow(temp_mother); hold on;
                % for x = 1: numel(mother_prop)
                %     plot(mother_prop(x).Centroid(1),mother_prop(x).Centroid(2),'ro');
                % end
                curr_mother_x = mother_prop(1).Centroid(1); %current mother mask's centroids
                curr_mother_y = mother_prop(1).Centroid(2);

               % (1) the previous y-centroid is 0, meaning there is no mother cell yet detected; and (2) the current y-centroid moves up within an acceptable range up and down from the previous y-centroid; and (3) the current x-centroid moves left/right within range
                y_up_allow = 10; %allowable vertical distance upwards from previous y-centroid
                y_down_allow = 75;
                x_allow = 18; %allow centroids to move horizontally up to 12px
                if (mother_y(i) == 0 | curr_mother_y < mother_y(i) + y_down_allow && curr_mother_y + y_up_allow >= mother_y(i) && abs(curr_mother_x-mother_x(i)) <= x_allow )
                    mother_x(i) = curr_mother_x; %update centroids
                    mother_y(i) = curr_mother_y;
                    mother_BW{i} = temp_mother; %update mother_BW for current col
                % if the current mother cell doesn't meet these criteria or there is no mother cell in the column, do not update the arrays and load the previous frame's mother_BW as the current temp_mother column mother mask
                else
                    temp_mother = mother_BW{i}; %restore previous frame's mother cell col mask without updating the centroids
                end
            else %if there are no mother cells/centroids, use previous mother mask
                temp_mother = mother_BW{i};
            end


            % % Add current column mask to the overall mask image for output
            % I_nuc_mask = horzcat(I_nuc_mask,BW3);
            % Add current mother cell mask to the overall mother mask image
            I_mother_mask = horzcat(I_mother_mask,temp_mother);


            % For each of the fluorescent channels...
            for y = 1:fluN
                curr_I_flu = I_flu{y}; % current fluorescent image
                I_flu_col = curr_I_flu(:,1+(i-1)*block:i*block); % current column in flu image

                % Determine the properties of the segmented cells using regionprops(BW_image,Intensity_Image,Properties)
                % Since we already have a column mask of just the single mother cell, we can simply track its PixelValues and there is no need to recalculate centroids
                col_prop = regionprops(temp_mother,I_flu_col,'PixelValues');

                % Structure for holding the cell fluorescence and other property data for all cells
                col_prop_2 = [col_prop(:).PixelValues];%the values of all the pixels in the mother cell mask

                % Top 80 fluorescence method #1: Take the top half of the values in the array
                col_prop_3 = sort(col_prop_2); %sort PixelValues in ascending order
                num_px = numel(col_prop_3); %number of pixels (area) of mother cell
                top_50 = floor(0.8*num_px+1):num_px; % top 50% range
                col_prop_4 = mean(col_prop_3(top_50)); %

                % % Top 50 fluorescence method #2: Find midpoint of min/max, and take the values above this threshold
                % mi_px = min(col_prop_2); ma_px = max(col_prop_2);
                % avg_px = (ma_px - mi_px)*0.5 + mi_px;
                % col_prop_3 = col_prop_2(col_prop_2 > avg_px);
                % col_prop_4 = mean(col_prop_3);

                % Store mother cell fluorescence in trajectories matrix
                traj(imid,i,y) = col_prop_4;
            end
        end

        % % Output mask image
        % imwrite(I_nuc_mask, out_name);
        % %figure; imshow(I_nuc_mask);

        % % Output mask+phase overlay
        % I_overlay = imfuse(I_ph,I_nuc_mask,'falsecolor','ColorChannels','red-cyan');
        % %figure; imshow(I_overlay)
        % imwrite(I_overlay, out_name_overlay)

        % Output mother mask+phase overlay
        I_overlay_2 = imfuse(I_ph,I_mother_mask,'falsecolor','ColorChannels','red-cyan');
        %figure; imshow(I_overlay_2);
        imwrite(I_overlay_2, out_name_mother)
    end

    % Write trajectory data to output_data
    save(output_data, 'traj');
end
