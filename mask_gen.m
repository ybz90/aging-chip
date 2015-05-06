function mask_gen(pos,imN)

    % mask_gen.m is used to generate a mask from the treshold and nuclear marker images. The threshold is used to determine the traps as columns, and within each one, the nuclear marker is segmented and dilated.
    % EdgeID is based on the MeanIntensity of the cell region.
    % (OPTIONAL) Export masks and mask+phase overlay images.
    % For the first frame: after identifying cells and cleaning their masks, the cell with the lowest y-position is identified as the mother cell. It's
    %


    % Original code by Meng Jin; last edited by Yuan Zhao 05/05/15

    % TO DO: INCLUDE CLEANING CODE FOR REMOVING NON CIRCULAR OBJECTS, ETC.
    % A potential fix regarding when the nuclear marker fades, resulting in a terrible segment vs background and failure of a mask; 1) do not allow for cells with a mask diameter <5 px?; 2) ignore potential segmented cells if no px > arbitrarily from manually searching 300au in intensity.


    % Create appropriate masks output directory (ie /masks/xy2526)
    mkdir(strcat('xy',pos,'/mask'))
    % Create output directory for mask overlaid atop phase images
    mkdir(strcat('xy',pos,'/mask_overlay'))



    outputdata =['output_xy',pos,'.mat'];

    % inialize empty 3D matrix to score trajectories (frame #, # flu channels, trap #)
    %Traj = zeros(imN, flrn, colN);




    for imid = 1:imN

        % input directory and image paths for stitched nuc and thrsh images
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        I_ph = imread(ph_name);

        thrsh_name = ['xy',pos,'/c1_thr/xy',pos,'_c1_thr_t',sprintf('%04g',imid),'.tif'];
        I0 = imread(thrsh_name);

        nuc_name = ['xy',pos,'/c3/xy',pos,'_c3_t',sprintf('%04g',imid),'.tif'];
        Iflr = imread(nuc_name);

        % output directory and mask image name
        save_name = ['xy',pos,'/mask/xy',pos,'_mask_t',sprintf('%04g',imid),'.tif'];
        save_name_c3ph = ['xy',pos,'/mask_overlay/xy',pos,'_mask_overlay_t',sprintf('%04g',imid),'.tif'];


        % initialize first frame
        if imid == 1

            % morphologically open threshold image, remove small objects less
            % than rect via rect struct element
            I1 = imopen(I0, strel('rectangle', [5, 100]));

            % possible that horizontal
            I1_extend = imdilate(I1, strel('line',400,0));
            I1_extend_lb = bwlabel(I1_extend);
            I1_lb = bwlabel(bwareaopen((I1_extend==0),5000));
            center_blck = (I1_lb==2);
            I1_cnt = (I1_lb==2).*(I0>0);

            % get vertical part
            I2 = imopen(I1_cnt, strel('rectangle',[10, 2]));

            % coarse grain vertical columns
            I3 = imdilate(I2, strel('rectangle',[100, 17]));
            I3_a = (I3>0) - imdilate((I1_extend_lb==1),strel('rectangle',[2,100]));
            I3_b = (I3_a>0);
            I3_c = bwareaopen(I3_b, 1500);
            I3_d = imdilate(I3_c, strel('disk',3));
            %figure; image(I3_c); %debug
            %figure; image(I3_d); %debug

            % label each trap object as a column [label matrix, num components]
            [I3_lb,colN]= bwlabeln(I3_c);

            Columns = I3_lb;
        end

        % initialize empty matrices of size ~ threshold image
        I_fulledge = zeros(size(I0));
        I_fillcell = I_fulledge;

        I_exm = I_fulledge;
        I_nuccell= I_fillcell;

        % init
        PAratio_nuc = cell(colN,1);

        fprintf('Generating mask. Frame number %d has %d columns.\n', imid, colN); %debug

        % for each column in the current frame
        for i=1:colN

            curr_col = (Columns ==i);
            Icf = double(Iflr).*curr_col; %double precision of nuc imread * current column

             % fluor prop of the column
            fplist_clm = regionprops(curr_col,Icf,'MeanIntensity','PixelIdxList');
            Iflr_clm_th = otsuthresh(Icf, fplist_clm.PixelIdxList); %otsu threshold current column nuc markers, calls otsuthresh.m

            %fill in holes between raw nuc image and otsu thresh; remove small objects less than 10px via area open
            Iflr_clm_mask0 = bwareaopen(imfill(Icf>Iflr_clm_th, 4, 'holes'),10);
            %dilate mask with disk
            Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));

            %label objects (cells)
            [Iflr_lb0, flr_N]= bwlabel(Iflr_mask_out);

            %if no objects (cells) are found, erode column and try again up to 5 times
            thrsh_drop=0;
            while flr_N==0 && thrsh_drop<5
                thrsh_drop = thrsh_drop+1;
                 temp_clm = curr_col;
                 column= imerode(temp_clm, strel('rectangle',[5,1]));
                 Icf = double(Iflr).*double(column);
                 fplist_clm = regionprops(curr_col,Icf,'MeanIntensity','PixelIdxList');
                Iflr_clm_th = otsuthresh(Icf, fplist_clm.PixelIdxList);
                Iflr_clm_mask0 = bwareaopen(imfill(Icf>(1-0.1*thrsh_drop)*Iflr_clm_th, 4, 'holes'),10);
                Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));

                [Iflr_lb0, flr_N]= bwlabel(Iflr_mask_out);

            end


            % Determine the properties of the labeled cells; regionprops(BW_image,Intensity_Image,Properties)
            prop_f1 = regionprops(Iflr_lb0,Icf,'Area','MeanIntensity','Centroid','MajorAxisLength', 'MinorAxisLength');

            % Get centroids of all cells and their x,y coordinates
            all_centroids = [prop_f1(:).Centroid];
            x_centroids = all_centroids(1:2:end-1);
            y_centroids = all_centroids(2:2:end);

            % Highest y-value of the centroids is the lowest cell in the trap (aka mother cell)
            % idx is its index in the centroids array; this is the ONLY cell we will track in propmx_f1 below
            [max_y idx] = max(y_centroids)

            % Structure for holding the cell fluorescence and other property data for all cells in the current column i
            % propmx_f1 has dimensions of 4 rows (props) x n cols (# of cells)
            % column idx stores the properties of the mother cell
            propmx_f1=[prop_f1(:).Area; prop_f1(:).MeanIntensity; x_centroids; y_centroids]

            fprintf('Column #%d has %d cells.\n', i, length(propmx_f1(:))/4); %debug


            % If at least one cell has been detected
            if ~isempty(propmx_f1)
                % two methods for connected cells:
                PAratio_nuc = propmx_f1(2,:).^2./propmx_f1(1,:)/4/pi;
                [nouse,clump_idx] = find(PAratio_nuc > 1.5);

                if ~isempty(clump_idx)
                    [PAratio_new, Iflr_mask2, Iflr_mask_label2] = de_clump_v2(Iflr_mask_out, Iflr_lb0, PAratio_nuc, clump_idx);
                    Iflr_mask_out = bwareaopen(Iflr_mask2, 10);
                end
            end

            I_nuccell = I_nuccell + Iflr_mask_out;


        end

        % Output mask image
        Iout = uint8((I_nuccell>0)*255);
        imwrite(Iout, save_name)
        %figure; imagesc(Iout) %debug

        % Output mask image overlaid on top of phase
        I_out2 = I_ph + uint16(I_nuccell)*2000;
        imwrite(I_out2, save_name_c3ph)
    end
    %clock

end
