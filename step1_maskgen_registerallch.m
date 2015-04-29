% Register all frames' positions in all 3 channels using the threshold image of the
% phase channel; the purpose of this code is to align all the frames to
% correct for drift and other movement over time

% Original code by Meng Jin; last edited by Yuan Zhao



% Input the position number to be registered; to align and account for
% camera movement and drift
pos = input('Input xy position: ');
pos = int2str(pos);

% Input the number of frames (times) to be processed
imN = input('Number of frames to process: ');

% In the directories for each channel and threshold, make a new subfolder 
% to store the registered files 
for channels = {'phase','thresh','flu','nuc'}
    mkdir(strcat('xy',pos,'/',channels{1},'/reg'))
end


% For every frame/timepoint
for imid = 0:imN-1 %assumes t starts at t0000
    
    % Load the names and output dirs of the raw frames to be registered,
    % for the current time frame (imid)
    thrsh_name = ['xy',pos,'/thresh/raw/xy',pos,'c1_pha_t',sprintf('%04g',imid),'.tif'];
    phase_name = ['xy',pos,'/phase/raw/xy',pos,'c1_t',sprintf('%04g',imid),'.tif'];
    fluor_name = ['xy',pos,'/flu/raw/xy',pos,'c2_t',sprintf('%04g',imid),'.tif'];
    nuclr_name = ['xy',pos,'/nuc/raw/xy',pos,'c3_t',sprintf('%04g',imid),'.tif'];
    
    % Set the names and output dirs of the registered frames for this imid
    thrsh_name2 = ['xy',pos,'/thresh/reg/xy',pos,'c1_pha_A_t',sprintf('%04g',imid),'.tif'];
    phase_name2 = ['xy',pos,'/phase/reg/xy',pos,'c1_A_t',sprintf('%04g',imid),'.tif'];
    fluor_name2 = ['xy',pos,'/flu/reg/xy',pos,'c2_A_t',sprintf('%04g',imid),'.tif'];
    nuclr_name2 = ['xy',pos,'/nuc/reg/xy',pos,'c3_A_t',sprintf('%04g',imid),'.tif'];
    
    % Import the frame images corresponding to the loaded names
    I0 = imread(thrsh_name); % Thresholding (I0) is used for alignments
    %figure, image(I0) %debug
    
    Iph = imread(phase_name); % Import phase, fluorescence, and nuclear so these images may also be registered
    Iflr = imread(fluor_name);
    Inuc = imread(nuclr_name);

    % Morphologically open binary threshold image w/ structuring element
    % (2x100 rect); ignore features in thresh < 2x75 rect
    % NOTE: Consider using gpuArray to accelerate this step?
    I1 = imopen(I0, strel('rectangle', [2, 75])); %changed from [2, 100]
    %figure, imshow(I1,[]); %debug
    I1_extend = imdilate(I1, strel('line',500,0)); %
    %figure, imshow(I1_extend,[]); %debug
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

    %imid %debug
    center_xy = regionprops(center_blck, 'Centroid');
    column_xy = regionprops(Ia_lb, 'Centroid'); 

    columnN = length(column_xy);
    curr_xy = round([column_xy(1).Centroid(1), center_xy.Centroid(2)]);

    % For first image, set as the current position
    if imid == 0
        use_refxy = curr_xy;

        Imask_new = I0;
        Iflr_new = Iflr;
        Inuc_new = Inuc;
        Iph_new  = Iph;
    
    % Update positions and shift  flu, nuc, and phase images
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
    imwrite(Imask_new, thrsh_name2);
    imwrite(Iflr_new, fluor_name2);
    imwrite(Inuc_new, nuclr_name2);
    imwrite(Iph_new, phase_name2);      
end

%clock