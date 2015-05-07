function mask_gen(pos,imN,colN)

    % mask_gen.m is used to generate a mask from the nuclear marker images.
    % The image is divided into colN slices for each trap as columns.
    % For each column: after identifying cells and cleaning their masks, the cell with the lowest y-position is identified as the mother cell.
    % Record the mother cell's position and fluorescence for all channels; and xy_pos and trap #.
    % (OPTIONAL) Export masks and mask+phase overlay images.

    % Yuan Zhao 05/06/2015


    % Create appropriate masks output directory (ie /masks/xy2526)
    mkdir(strcat('xy',pos,'/mask'))
    % Create output directory for mask overlaid atop phase images
    mkdir(strcat('xy',pos,'/mask_overlay'))



    % inialize empty 3D matrix to score trajectories (frame #, # flu channels, trap #)
    %Traj = zeros(imN, flrn, colN);

    % Output trajectory data for all traps on this xy position
    outputdata =['output_xy',pos,'.mat'];



    % define column slice width as the image width divided by  colN
    block = round(512/colN);

    for imid = 1:imN

        % input directory and image paths for nuclear marker, fluorescence, and phase channel images
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph = imread(ph_name);

        nuc_name = ['xy',pos,'/c3/xy',pos,'_c3_t',sprintf('%04g',imid),'.tif'];
        Inuc = imread(nuc_name);

        % output directory and mask image name
        save_name = ['xy',pos,'/mask/xy',pos,'_mask_t',sprintf('%04g',imid),'.tif'];
        save_name_c3ph = ['xy',pos,'/mask_overlay/xy',pos,'_mask_overlay_t',sprintf('%04g',imid),'.tif'];

        %fprintf('Generating mask. Frame number %d has %d columns.\n', imid, colN); %debug

        % initialize empty image for storing mask
        Inuc_mask = [];

        % for each column in the current frame
        for i = 1:colN

            Icf = Inuc(:,1+(i-1)*block:i*block);
            %figure; imagesc(Inuc(:,1+(i-1)*block:i*block)) %debug

            [level EM] = graythresh(Icf);
            BW = im2bw(Icf,0.05);


            %fill in holes between raw nuc image and otsu thresh; remove small objects less than 10px via area open
            BW2 = imfill(BW,'holes');
            %dilate mask with disk
            BW3 = imdilate(BW2, strel('disk',1));

            %figure; imshow(BW3) %debug
            Inuc_mask = horzcat(Inuc_mask,BW3);

            % CODE FOR DECLUMPING MAYBE FOR CONNETED CELLS? IDK IF NEEDED






            % Determine the properties of the labeled cells; regionprops(BW_image,Intensity_Image,Properties)
            prop_f1 = regionprops(BW3,Icf,'Area','MeanIntensity','Centroid','MajorAxisLength', 'MinorAxisLength');


            % Get centroids of all cells and their x,y coordinates
            all_centroids = [prop_f1(:).Centroid];
            x_centroids = all_centroids(1:2:end-1);
            y_centroids = all_centroids(2:2:end);

            % Structure for holding the cell fluorescence and other property data for all cells in the current column i
            % propmx_f1 has dimensions of 4 rows (props) x n cols (# of cells)
            % column idx stores the properties of the mother cell
            propmx_f1=[prop_f1(:).Area; prop_f1(:).MeanIntensity; x_centroids; y_centroids]
            %fprintf('Column #%d has %d cells.\n', i, length(propmx_f1(:))/4); %debug

            % Highest y-value of the centroids is the lowest cell in the trap (aka mother cell)
            % idx is the column # of the cell we want to track in propmx_f1
            [max_y idx] = max(y_centroids);


            % TRAJ STUFF


        end

        %figure; imshow(Inuc_mask)

        % Output mask image
        %Iout = uint8((I_nuccell>0)*255);
        %imwrite(Iout, save_name)
        %figure; imagesc(Iout) %debug

        % Output mask image overlaid on top of phase
        %I_out2 = I_ph + uint16(I_nuccell)*2000;
        %imwrite(I_out2, save_name_c3ph)
    end
    %clock

end
