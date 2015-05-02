% Mask generation code to generate masks from registered images

% based on v6; in v6, edgeid is based on MeanIntensity of cell region.
% in v7, considering glu 0.02 when cells are more ellipsoid like shapes.
% and connected

% Original code by Meng Jin; last edited by Yuan Zhao


clear all
%clock

% Input the stitched position name (ie 2526)
pos = int2str(input('Input xy position: '));

% Input the number of time frames for which to generate masks
imN = input('Number of time frames: ');

% Create appropriate masks output directory (ie /masks/xy2526)
mkdir(strcat('xy',pos,'/mask_raw'))


for imid = 0:imN-1
    
    % imid %debug
    
    % input directory and image paths for stitched nuc and thrsh images
    phase_name = ['xy',pos,'/c1_thr/reg/xy',pos,'_c1_thr_reg_t',sprintf('%04g',imid),'.tif'];
    nuc_name = ['xy',pos,'/c3/reg/xy',pos,'_c3_reg_t',sprintf('%04g',imid),'.tif'];
    
    % output directory and mask image name
    save_name = ['xy',pos,'/mask_raw/xy',pos,'_mask_raw_t',sprintf('%04g',imid),'.tif'];

    I0 = imread(phase_name);
    Iflr = imread(nuc_name);
    

    if imid == 0

        I1 = imopen(I0, strel('rectangle', [3, 100])); %changed from [5, 100]

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

        [I3_lb colN]= bwlabeln(I3_c);

        Columns = I3_lb;
    end
    
    I_fulledge = zeros(size(I0));
    I_fillcell = I_fulledge;

    I_exm = I_fulledge;
    I_nuccell= I_fillcell;

    PAratio_nuc = cell(colN,1);

        if colN ~= 21;
            fprintf('image %d has %d columns.\n', imid, colN);
        end

    for i=1:colN

        curr_col = (Columns ==i);
        Icf = double(Iflr).*curr_col;

         % fluor prop of the column
        fplist_clm = regionprops(curr_col,Icf,'MeanIntensity','PixelIdxList'); 
        Iflr_clm_th = otsuthresh(Icf, fplist_clm.PixelIdxList);
        Iflr_clm_mask0 = bwareaopen(imfill(Icf>Iflr_clm_th, 4, 'holes'),10);
        Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));

        [Iflr_lb0, flr_N]= bwlabel(Iflr_mask_out);

        thrsh_drop=0;
        while flr_N==0
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

        prop_f1 = regionprops(Iflr_lb0,Icf,'Area','MeanIntensity','Solidity','Perimeter');


        if max([prop_f1.Area])>500
            Iflr_clm_mask0 = bwareaopen(imfill(Icf>1.1*Iflr_clm_th, 4, 'holes'),10);
            Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));

            Iflr_lb0= bwlabel(Iflr_mask_out);
            prop_f1 = regionprops(Iflr_lb0,Icf,'Area','MeanIntensity','Solidity','Perimeter');
        end

        propmx_f1=[prop_f1(:).Area ; prop_f1(:).Perimeter; prop_f1(:).Solidity; prop_f1(:).MeanIntensity];


        if ~isempty(propmx_f1)
            % two methods for connected cells:
            PAratio_nuc = propmx_f1(2,:).^2./propmx_f1(1,:)/4/pi;
            [nouse clump_idx] = find(PAratio_nuc > 1.5); 

            if ~isempty(clump_idx)
                [PAratio_new, Iflr_mask2, Iflr_mask_label2] = de_clump_v2(Iflr_mask_out, Iflr_lb0, PAratio_nuc, clump_idx);
                Iflr_mask_out = bwareaopen(Iflr_mask2, 10);
            end
        end

        I_nuccell = I_nuccell + Iflr_mask_out;
    end
    
% figure; imagesc(I_nuccell)

    Iout = uint8((I_nuccell>0)*255);
    imwrite(Iout, save_name)
end

%clock