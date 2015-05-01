% Mask generation code to generate masks from stitched nuclear marker and
% threshold images for cell tracking

% Original code by Meng Jin; last edited by Yuan Zhao

% based on v6; in v6, edgeid is based on MeanIntensity of cell region.
% in v7, considering glu 0.02 when cells are more ellipsoid like shapes.
% and connected



clear all
clock

% Input the stitched position name (ie 2526)
mergexy = int2str(input('Input composite xy name: '));

% Input the number of time frames for which to generate masks
imN = input('Number of time frames: ');

% Create appropriate masks output directory (ie /masks/xy2526)
mkdir(strcat('xy',mergexy,'/mask_raw'))


for imid = 0:imN-1
    
    imid %debug
    
    % input directory and image paths for stitched nuc and thrsh images
    phase_name = ['xy',mergexy,'/thresh/xy',mergexy,'c1_pha_t',sprintf('%04g',imid),'.tif'];
    nuc_name = ['xy',mergexy,'/nuc/xy',mergexy,'c3_t',sprintf('%04g',imid),'.tif'];
    
    % output directory and mask image name
    save_name = ['xy',mergexy,'/mask_raw/ph_mask_raw_t',sprintf('%04g',imid),'.tif'];

    
I0 = imread(phase_name);
Iflr = imread(nuc_name);
% Inuc0 = imread('test_nuc_a0125.tif');

I1 = imopen(I0, strel('rectangle', [2, 100])); %changed from 5, 100
% figure; imagesc(I1)

% possible that horizontal 
I1_extend = imdilate(I1, strel('line',400,0)); %changed from 400,0
I1_extend_lb = bwlabel(I1_extend);

% imagesc(I1_extend)

I1_lb = bwlabel(bwareaopen((I1_extend==0),5000));
center_blck = (I1_lb==2);

 
I1_cnt = (I1_lb==2).*(I0>0);


% get vertical part
I2 = imopen(I1_cnt, strel('rectangle',[10, 2]));

% coarse grain vertical columns
I3 = imdilate(I2, strel('rectangle',[100, 17]));
I3_a = (I3>0) - imdilate((I1_extend_lb==1),strel('rectangle',[2,100])); %changed from 2, 100
I3_b = (I3_a>0);
I3_c = bwareaopen(I3_b, 1500);

[I3_lb colN ]= bwlabeln(I3_c);

Ia = I3.*I1_cnt;


Ia1 = imdilate(Ia, strel('disk',6)); %
Ia2 = imfill(Ia1, 'holes');
Ia3 = bwareaopen(Ia2, 1000); % remove any non-column features

% [Ia_lb colN ]= bwlabeln(Ia3);


I_fulledge = zeros(size(I0));
I_fillcell = I_fulledge;

I_exm = I_fulledge;
% use input box/column to get local column


I_fulledge = zeros(size(I0));
I_fillcell = I_fulledge;

I_nuccell= I_fillcell;

PAratio_nuc = cell(colN,1);
for i=1:colN
%     Ic1 = double(I0).*double(Ia_lb ==i); 
%     Ic1b = bwareaopen(Ic1, 25);
    Icf = double(Iflr).*double(I3_lb==i);
    
%     temp1 = imerode((Ia_lb ==i),strel('disk',4)); % get a thinner column block (from imdilate(disk,5))
%     temp2 = imerode(temp1,strel('disk',4));
%     tempe = temp1-temp2; % get the edge of each column
    
%     Ic2 = tempe + (Ic1b>0);
%     Ic3 = bwareaopen(Ic2, 20);
%     
%     I_fulledge = I_fulledge + Ic3;
%     I_cell = bwareaopen((Ia_lb ==i) - Ic3, 10); 
%     [I_cell_lb, celln] = bwlabel(I_cell);
     
     
     % fluor prop of the column
    fplist_clm = regionprops((I3_lb==i),Icf,'MeanIntensity','PixelIdxList'); 
    Iflr_clm_th = otsuthresh(Icf, fplist_clm.PixelIdxList); %calls otsu's thresholding algorithm function otsuthresh.m
    Iflr_clm_mask0 = bwareaopen(imfill(Icf>Iflr_clm_th, 4, 'holes'),10);
    Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));
    
    [Iflr_lb0, flr_N]= bwlabel(Iflr_mask_out);
    
    while flr_N==0
         temp_clm = (I3_lb==i);
         column= imerode(temp_clm, strel('rectangle',[5,1]));
         Icf = double(Iflr).*double(column);
         fplist_clm = regionprops((I3_lb==i),Icf,'MeanIntensity','PixelIdxList'); 
        Iflr_clm_th = otsuthresh(Icf, fplist_clm.PixelIdxList);
        Iflr_clm_mask0 = bwareaopen(imfill(Icf>0.9*Iflr_clm_th, 4, 'holes'),10);
        Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));

        [Iflr_lb0, flr_N]= bwlabel(Iflr_mask_out);
        
    end
    
    prop_f1 = regionprops(Iflr_lb0,Icf,'Area','MeanIntensity','Solidity','Perimeter');
    
    
    if max([prop_f1.Area])>500
        Iflr_clm_mask0 = bwareaopen(imfill(Icf>1.25*Iflr_clm_th, 4, 'holes'),10);
        Iflr_mask_out = imdilate(Iflr_clm_mask0, strel('disk',1));
    
        Iflr_lb0= bwlabel(Iflr_mask_out);
        prop_f1 = regionprops(Iflr_lb0,Icf,'Area','MeanIntensity','Solidity','Perimeter');
    end
    
    propmx_f1=[prop_f1(:).Area ; prop_f1(:).Perimeter; prop_f1(:).Solidity; prop_f1(:).MeanIntensity];
    
    % two methods for connected cells:
    PAratio_nuc = propmx_f1(2,:).^2./propmx_f1(1,:)/4/pi;
    
    [nouse clump_idx] = find(PAratio_nuc > 1.5); 
    
     
    if ~isempty(clump_idx)
        [PAratio_new, Iflr_mask2, Iflr_mask_label2] = de_clump_v2(Iflr_mask_out, Iflr_lb0, PAratio_nuc, clump_idx);

        Iflr_mask_out = bwareaopen(Iflr_mask2, 10);
    end
    
    I_nuccell = I_nuccell + Iflr_mask_out;
end
 
% figure; imagesc(I_nuccell)
    
    Iout = uint8((I_nuccell>0)*255);
    imwrite(Iout, save_name)
     
end

clock