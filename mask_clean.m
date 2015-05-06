function mask_clean(pos,imN)

    % mask_clean.m is used to further clean the previously generated masks, using a cytosolic fluorescent channel to remove dead cells and non-cell debris/artifacts. This function will export the further cleaned masks and export an image of it overlayed on top of the phase images.

    % Original code by Meng Jin; last edited by Yuan Zhao 05/05/15


    % Create output directories for cleaned mask images, and also overlayed
    % phase images
    mkdir(strcat('xy',pos,'/mask_clean'))
    mkdir(strcat('xy',pos,'/mask_clean_phase'))


    for imid = 1:imN

        imid %debug

        % input image names
        ph_name = ['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif'];
        thrsh_name = ['xy',pos,'/c1_thr/xy',pos,'_c1_thr_t',sprintf('%04g',imid),'.tif'];

        fluor_name = ['xy',pos,'/c2/xy',pos,'_c2_t',sprintf('%04g',imid),'.tif'];

        mask_name = ['xy',pos,'/mask_raw/xy',pos,'_mask_raw_t',sprintf('%04g',imid),'.tif'];

        % output image names
        save_name_mask = ['xy',pos,'/mask_clean/xy',pos,'ph_mask_clean_t',sprintf('%04g',imid),'.tif'];
        save_name_c3ph = ['xy',pos,'/mask_clean_phase/xy',pos,'ph_mask_clean_phase_t',sprintf('%04g',imid),'.tif'];


        I0= imread(thrsh_name);
        Iflr = imread(fluor_name);
        Imask = imread(mask_name);
        I_ph = imread(ph_name);

        % initialize the column mask
        if imid == 0

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

            % label connected components in image; vertical columns [Label
            % matrix, num components]
            [I3_lb,colN] = bwlabeln(I3_c);

            Columns = I3_lb;
        end


        for i = 1:colN

            curr_col = (Columns ==i);

            Ic2 = double(Iflr).*curr_col;  % wholecell fluor in this column
            Ic3 = double(Imask).*curr_col; % nuclear mask
            [mask_lb, nucN] = bwlabel(Ic3);

            if nucN>0
                col_c2 = regionprops(curr_col,Ic2,'MeanIntensity','PixelIdxList');
                col_c3 = regionprops(mask_lb,Ic2,'MeanIntensity','PixelIdxList');

                cutoff_c2 = otsuthresh(Ic2,col_c2.PixelIdxList);

                Ic2_inmask = [col_c3(:).MeanIntensity];

                remove_maskid = find(Ic2_inmask < cutoff_c2);

                if ~isempty(remove_maskid)
                    for kk = 1:length(remove_maskid)
                        remove_id = remove_maskid(kk);
                        temp_mask = (mask_lb == remove_id);
                        new_mask = mask_lb - remove_id*temp_mask;
                        mask_lb = new_mask;
                    end

                    Imask = (Imask>0) - (Ic3>0) + (mask_lb>0);
                end
            end
        end

        Imask2 = bwareaopen((Imask>0), 14);


        Iout = uint8((Imask2>0)*255);
        imwrite(Iout, save_name_mask)


        I_out2 = I_ph + uint16(Imask2)*2000;
        imwrite(I_out2, save_name_c3ph)
    end


end
