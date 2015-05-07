function mask_traj(pos,imN,colN,fluN)

    % mask_traj.m is used to generate a mask from the nuclear marker images and produce trajectories for the mother cells in each trap.
    % The image is divided into colN slices for each trap as columns.
    % Each column has a mask generated, and the cell with the lowest y-coord is identified as the mother cell.
    % Record the mother cell's position and fluorescence for all channels; and xy_pos and trap #.
    % Export masks and mask+phase overlay images.

    % Yuan Zhao 05/06/2015


    % Create output directory for storing masks
    mkdir(strcat('xy',pos,'/mask'))
    % Create output directory for storing mask+phase overlays
    mkdir(strcat('xy',pos,'/mask_overlay'))

    % Initialize empty 3D matrix to store trajectories (frame #, trap #, flu channel #)
    traj = zeros(imN, colN, fluN);
    % Output trajectory data for all traps on this xy position
    outputdata =['xy',pos,'/output_xy',pos,'.mat']; %include position and date in the filename

    % Define column slice width as the image width divided by colN
    block = round(512/colN);

    % For each frame...
    for imid = 1:imN
        %fprintf('Generating mask. Frame number %d has %d columns.\n', imid, colN); %debug

        % Input directory and image paths for nuclear marker, fluorescence, and phase channel images
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph = imread(ph_name);

        nuc_name = ['xy',pos,'/c3/xy',pos,'_c3_t',sprintf('%04g',imid),'.tif'];
        I_nuc = imread(nuc_name);


        % IMPORT FLUORESENCE CHANNEL(S) FOR MEASURING


        % Output directory and mask image name
        out_name = ['xy',pos,'/mask/xy',pos,'_mask_t',sprintf('%04g',imid),'.tif'];
        out_name_overlay = ['xy',pos,'/mask_overlay/xy',pos,'_mask_overlay_t',sprintf('%04g',imid),'.tif'];

        % Initialize mask image for this frame
        %I_nuc_mask = [];
        I_nuc_mask = zeros(512,1); %since 512 is not divisible by 7, the block size (col width) is rounded and we need to add a column of black pixels to get 511+1 = 512 px

        % For each column in the current frame...
        for i = 1:colN


            % ADD CODE TO CHECK IF THERE ANY CELLS AT ALL; IGNORE TRAP IF NONE FOUND


            Icf = I_nuc(:,1+(i-1)*block:i*block); %current column of nuclear marker image
            %figure; imagesc(I_nuc(:,1+(i-1)*block:i*block)) %debug

            % Otsu's threshold nuclear marker image to obtain binary cell image mask
            [level EM] = graythresh(Icf);
            BW = im2bw(Icf,0.05);
            BW2 = imfill(BW,'holes');
            BW3 = imdilate(BW2, strel('disk',1)); %dilate mask with disk
            %figure; imshow(BW3) %debug


            % IMPLEMENT CODE FOR DECUMPLING CONNECTED CELLS AS WELL AS REMOVING ON CIRCULAR CELLS; this is unimportant if the lower-most cell is unaffected by these kind of mask segmentation issues


            I_nuc_mask = horzcat(I_nuc_mask,BW3); %add current column mask to the overall mask image for output

            % Determine the properties of the segmented cells using regionprops(BW_image,Intensity_Image,Properties)
            col_prop = regionprops(BW3,Icf,'Area','MeanIntensity','Centroid','MajorAxisLength', 'MinorAxisLength');

            % Get centroids of all cells and their x,y coordinates
            all_centroids = [col_prop(:).Centroid];
            x_centroids = all_centroids(1:2:end-1);
            y_centroids = all_centroids(2:2:end);

            % Structure for holding the cell fluorescence and other property data for all cells
            % col_prop_2 has dimensions of 4 rows (cell properties) x n cols (# of cells in the trap)
            col_prop_2 = [col_prop(:).Area; col_prop(:).MeanIntensity; x_centroids; y_centroids];
            %fprintf('Column #%d has %d cells.\n', i, length(col_prop_2(:))/4); %debug

            % The centroid with the highest y-coordinate value is the cell that is lowest in the trap, aka the mother cell; its index, idx, specifies the column in col_prop_2 that stores its properties
            [max_y idx] = max(y_centroids);

            % Properties of the mother cell, from the column #idx
            mother_prop = col_prop_2(:,idx);

            % Store mother cell fluorescence in trajectories matrix
            traj(imid,i,1) = mother_prop(2)
        end

        % Output mask image
        imwrite(I_nuc_mask, out_name)
        figure; imshow(I_nuc_mask)

        % Output mask+phase overlay
        C = imfuse(I_nuc_mask,I_ph,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
        imwrite(C, out_name_overlay)
        %figure; imshow(C)
    end

    % WRITE TRAJECTORY TO OUTPUT .MAT FILE
end
