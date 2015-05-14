function mask_traj(pos,imN,colN,fluN)

    % mask_traj.m is used to generate a mask from the nuclear marker images and produce trajectories for the mother cells in each trap.
    % The image is divided into colN slices for each trap as columns.
    % Each column has a mask generated, and the cell with the lowest y-coord is identified as the mother cell.
    % Record the mother cell's position and fluorescence for all channels; and xy_pos and trap #.
    % Export masks and mask+phase overlay images.

    % Yuan Zhao 05/13/2015


    % Mask generation

    % Create output directory for storing masks
    mkdir(strcat('xy',pos,'/mask'))
    % Create output directory for storing mask+phase overlays
    mkdir(strcat('xy',pos,'/mask_overlay'))
    % Create output directory for storing mother only mask+phase overlays
    mkdir(strcat('xy',pos,'/mother_mask_overlay'))

    % Get image dimensions
    ph_1 = ['xy',pos,'/c1/xy',pos,'_c1_t0001.tif'];
    I_ph_1 =imread(ph_1);
    I_size = size(I_ph_1); % rows (image height), cols (image width)
    height = I_size(1);
    width = I_size(2);

    % Define column slice width as the image width divided by colN
    block = round(width/colN);

    % Initialize overall overlaid mask of mother cell across all frames
    overall_mother = zeros(height,width); %%%%%%%%%%%%%%
    overall_daughter = zeros(height,width);

    % Store the centroids of each mother cell during early frames
    mother_x = cell(1,colN);
    mother_y = cell(1,colN);


    % For each frame...
    for imid = 1:imN

        fprintf('Generating mask for frame number %d.\n', imid); %debug

        % Input directory and image paths for nuclear marker and phase channel images
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph = imread(ph_name);

        nuc_name = ['xy',pos,'/c3/xy',pos,'_c3_t',sprintf('%04g',imid),'.tif'];
        I_nuc = imread(nuc_name);

        % Output directory and mask image names
        out_name = ['xy',pos,'/mask/xy',pos,'_mask_t',sprintf('%04g',imid),'.tif'];
        out_name_overlay = ['xy',pos,'/mask_overlay/xy',pos,'_mask_overlay_t',sprintf('%04g',imid),'.tif'];


        % Initialize mask images for this frame
        I_size = size(I_nuc);
        if rem(width,colN) == 0 % if image width is divisible by 7, init blank mask array
            I_nuc_mask = []; % mask of all cells
            I_mother_mask = []; % mask of mother cell only
            I_daughter_mask = []; % mask of lowest daughter cell only
        else
            I_nuc_mask = zeros(height,1); % since 512 is not divisible by 7, the block size (col width) is rounded and we need to add a column of black pixels to get 511+1 = 512 px
            I_mother_mask = zeros(height,1); %
            I_daughter_mask = zeros(height,1); %
        end


        % For each column in the current frame...
        for i = 1:colN

            Icf = I_nuc(:,1+(i-1)*block:i*block); %current column of nuclear marker image
            %figure; imagesc(Icf) %debug

            % Otsu's threshold nuclear marker image to obtain binary column mask
            %[level EM] = graythresh(Icf);
            level = 0.05;
            BW = im2bw(Icf,level);
            BW2 = imfill(BW,'holes');
            BW3 = imdilate(BW2, strel('disk',1)); %dilate mask with disk
            %figure; imshow(BW3) %debug





            % INSERT mask_optimize.m CODE FRAGMENT HERE

            % I have commented out this section of mask optimization in the interest of speed and performance. Since all of the masks are simply being overlaid atop one another anyway, it is a reasonable presumption that there will be little difference if a few masks are missing on certain frames, versus if each frame's mask is used individually, in which case it is crucial that the mother cell is correctly detected and thresholded.


            mask_prop = regionprops(BW3,Icf,'Area','Centroid'); %areas and centroids of the cell masks for the current column

            % Get centroids of all cells and their x,y coordinates
            all_centroids = [mask_prop(:).Centroid];
            x_centroids = all_centroids(1:2:end-1);
            y_centroids = all_centroids(2:2:end);

            % Structure for holding the cell area and centroid data for all cells
            % mask_prop_2 has dimensions of 3 rows (cell properties) x n cols (# of cells in the trap)
            mask_prop_2 = [mask_prop(:).Area; x_centroids; y_centroids];
            %fprintf('Column #%d has %d cells.\n', i, length(mask_prop_2(:))/4); %debug

            % The centroid with the highest y-coordinate value is the cell that is lowest in the trap, aka the mother cell; its index, idx, specifies the column in mask_prop_2 that stores its properties
            [max_y,idx] = max(y_centroids);

            % Properties of the mother cell, from the column #idx
            mother_prop = mask_prop_2(:,idx); %mother_prop(1) is the area of the mother cell
            areas = mask_prop_2(1,:);


            % Mask optimization code for if cells are not detected or if they do not meet certain quality conditions
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
            %NOTE: Remove all cells below a certain area, to remove the artifacts.
            %NOTE: Add code for declumping cells; since we only care about the lowest mother cell, this will only be needed for cases where a mother and a daughter are still attached when the frame was taken, and so we can approach this using a noncircularity method to identify these scenarios. Then we can do binary watershed or something similar and take the lower most object as the mother.





            % bwlabel does column-wise search by default; so to do row-wise searching for the lowest object, we transpose the BW3 binary image input and then transpose back the output of bwlabel
            BW4 = BW3.';
            [L,num] = bwlabel(BW4); %

            if num == 0 %if there is no cells in the column, use a pure black, empty column (all zeros)
                temp_mother = zeros(height,block);
            else
                temp_mother = (L==num); %image is only the mask of the mother cell
                temp_mother = temp_mother.';
            end

            if num-1 == 0 %if there are no daughter cells in the column, use a pure black, empty column (all zeros)
                temp_daughter = zeros(height,block);
            else
                temp_daughter = (L==num-1);
                temp_daughter = temp_daughter.';
            end


            % Add current mother and nearest daughter cell masks to this frame's overall mask images
            I_mother_mask = horzcat(I_mother_mask,temp_mother);
            I_daughter_mask = horzcat(I_daughter_mask,temp_daughter);

            % Add current column mask to the overall mask image for output
            I_nuc_mask = horzcat(I_nuc_mask,BW3);


            % Find centroid of mother cell, add to mother_x and mother_y arrays
            num_train = 100; % the number of frames during which to train the initial centroid coordinates
            if imid <= num_train
                mother_prop = regionprops(temp_mother,'centroid'); % find centroid of the current column mother cell
                if ~isempty(mother_prop) %if there is a mother cell in this particular column
                    %numel(mother_prop)
                    mother_prop2 = [mother_prop(1).Centroid(1),mother_prop(1).Centroid(2)];
                    % figure; imshow(temp_mother); hold on;
                    % for x = 1: numel(mother_prop)
                    %     plot(mother_prop(x).Centroid(1),mother_prop(x).Centroid(2),'ro');
                    % end
                    mother_x{i}(imid) = mother_prop(1).Centroid(1);
                    mother_y{i}(imid) = mother_prop(1).Centroid(2);
                end
            end
        end

        % Output current frame mask image
        imwrite(I_nuc_mask, out_name)
        %figure; imshow(I_nuc_mask)

        % Output current frame mask+phase overlay
        I_overlay = imfuse(I_ph,I_nuc_mask,'falsecolor','ColorChannels','red-cyan');
        %figure; imshow(I_overlay)
        imwrite(I_overlay, out_name_overlay)

        % Overlay current frame mother, daughter cell masks on top of previous merged overlay
        overall_mother = overall_mother + I_mother_mask;
        overall_daughter = overall_daughter + I_daughter_mask;
    end

    mother_only_mask = overall_mother - overall_daughter;
    %figure; imshow(overall_mother)
    %figure; imshow(overall_daughter)
    mother_only_mask = mother_only_mask>0;
    mother_only_mask = bwareaopen(mother_only_mask,25); %remove objects smaller than 25px

    % To clean the mask, we will remove all mask pixels in the binary image outside of a square boundary of length = 20(?) px enclosing the average centroid of each mother cell's first 50(?) frames.
        % To do this, we will 2*mother_only_mask, to make 1 values into 2. Next, we will subtract the binary image square from the 2*mask, resulting in values of:
        % 2 for mask outside of square; 1 for mask overlapping with square; -1 for square outside of box, and 0 for empty space
        % Taking only pixels with value 1 gives us the mask that overlaps with the square

    % Generate square binary image with which to constrict mother_only_mask
    constrictBW = zeros(height,width);
    for k = 1:colN % for every centroid
        x_cen = mean(mother_x{k})+(k-1)*block; %add block*(k-1) to x-coord account for column of centroid
        x_cen = ceil(x_cen);
        y_cen = mean(mother_y{k});
        y_cen = ceil(y_cen);
        s = 8; % length of square's side
        for xc = x_cen-s:x_cen+s % create a square around centroid of side length 2*s
            for yc = y_cen-s:y_cen+s
                %[xc,yc];
                constrictBW(yc,xc) = 1;
            end
        end
    end
    %figure; imshow(constrictBW);
    %figure; imshow(mother_only_mask);

    % Subtract constrictBW from 2*mother_only_mask
    double_mother = 2*mother_only_mask;
    mother_constrict = double_mother - constrictBW;

    % Generate final mother cell mask, constrained by centroid-centered squares
    mother_only_final = zeros(height,width);
    for l = 1:height
        for m = 1:width
            if mother_constrict(l,m) == 1
                mother_only_final(l,m) = 1;
            end
        end
    end
    %figure; imshow(mother_only_final);

    % Save mother_only_final combined mask
    imwrite(mother_only_final, ['xy',pos,'/xy',pos,'_overall_mask.tif']);


    % Trajectory calculation

    % Initialize empty 3D matrix to store trajectories (frame #, trap #, flu channel #)
    traj = zeros(imN, colN, fluN);
    % Output trajectory data for all traps on this xy position
    output_data =['xy',pos,'/xy',pos,'_traj.mat']; %include position and date in the filename

    % Label all the mother cells in the mother_only_final, from left to right (col-wise is default)
    [M,N] = bwlabel(mother_only_final);
    %figure; imshow(mother_only_final);

    % For each frame...
    for imid = 1:imN

        fprintf('Calculating trajectory for frame number %d.\n', imid); %debug

        % Input directory and image paths for every fluorescent channel
        flu_names = cell(1,fluN);
        I_flu = cell(1,fluN);
        for z = 1:fluN
            flu_names{z} = ['xy',pos,'/c',num2str(1+z),'/xy',pos,'_c',num2str(1+z),'_t',sprintf('%04g',imid),'.tif'];
            I_flu{z} = imread(flu_names{z});
        end

        % For each mother cell / column...
        for j = 1:N
            curr_mother = (M==j);
            %figure; imshow(curr_mother)

            % For each of the fluorescent channels...
            for y = 1:fluN
                curr_I_flu = I_flu{y};
                %figure; imshow(curr_I_flu)

                % Determine the properties of the mother cell using regionprops(BW_image,Intensity_Image,Properties)
                col_prop = regionprops(curr_mother,curr_I_flu,'PixelValues');

                % Structure for holding the cell fluorescence and other property data
                col_prop_2 = [col_prop(:).PixelValues]; %the values of all the pixels in the mother cell mask
                % col_prop_3 = sort(col_prop_2); %sort PixelValues in ascending order
                % num_px = numel(col_prop_3); %number of pixels (area) of mother cell
                % top_50 = floor(0.5*num_px+1):num_px; % top 50% range
                % col_prop_4 = mean(col_prop_3(top_50)); %
                mi_px = min(col_prop_2); ma_px = max(col_prop_2);
                avg_px = (ma_px - mi_px)*0.5 + mi_px;
                col_prop_3 = col_prop_2(col_prop_2 > avg_px);
                col_prop_4 = mean(col_prop_3);

                % Store mother cell fluorescence in trajectories matrix
                traj(imid,j,y) = col_prop_4;
            end
        end
    end


    % Export mother_only_final overlayed on top of phase
    for imid = 1:imN
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph = imread(ph_name);

        out_name_mother_overlay = ['xy',pos,'/mother_mask_overlay/xy',pos,'_mother_mask_overlay_t',sprintf('%04g',imid),'.tif'];

        I_mother_ph = imfuse(I_ph,mother_only_final,'falsecolor','ColorChannels','red-cyan');
        imwrite(I_mother_ph, out_name_mother_overlay)
    end

    % Write trajectory data to output_data
    save(output_data, 'traj');
end
