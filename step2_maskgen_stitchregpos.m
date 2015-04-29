% This code stitches (horizontally concatenates) xy positions together to
% create a single image for each time frame; to be used on output images
% after registering them to correct for positional drift

% Original code by Meng Jin; last edited by Yuan Zhao



% Input xy position #s of images to be stitched
pos_input = input('Positions to stitch {1, ... n}: ');

% Convert input positions from ints to strs
pos = {};
for k = 1:length(pos_input)
    pos{k} = int2str(pos_input{k});
end
posn = length(pos);

% Input desired output name of stitched image
outname = int2str(input('Output composite xy name: '));

% Make output directories named 'stitched/xy_outname/'
for channels = {'phase','thresh','flu','nuc'}
    mkdir(strcat('xy',outname,'/',channels{1}))
end

% Input the number of time frames to stitch
imN = input('Number of time frames: ');

%create 1xn array for storing filenames
thrsh_name = cell(1, posn+1); 
phase_name = cell(1, posn+1);
fluor_name = cell(1, posn+1);
nuclr_name = cell(1, posn+1);

for imid = 0:imN-1 %for each time frame, assumes t starts at t0000
    
    %store file names for each of the positions to stitch in array indices 1 to n-1
    %for binary threshold, phase, and 2 fluorescent channels
    for j = 1:posn
        thrsh_name{j} = ['xy',pos{j},'/thresh/reg/xy',pos{j},'c1_pha_A_t',sprintf('%04g',imid),'.tif'];
        phase_name{j} = ['xy',pos{j},'/phase/reg/xy',pos{j},'c1_A_t',sprintf('%04g',imid),'.tif'];
        fluor_name{j} = ['xy',pos{j},'/flu/reg/xy',pos{j},'c2_A_t',sprintf('%04g',imid),'.tif'];
        nuclr_name{j} = ['xy',pos{j},'/nuc/reg/xy',pos{j},'c3_A_t',sprintf('%04g',imid),'.tif'];
    
        %thrsh_name{j} =['threshold/reg/xy',pos{j},'/xy',pos{j},'c1_pha_A_t',sprintf('%04g',imid),'.tif'];
        %phase_name{j} =['phase/reg/xy',pos{j},'/xy',pos{j},'c1_A_t',sprintf('%04g',imid),'.tif'];
        %fluor_name{j} =['flu/reg/xy',pos{j},'/xy',pos{j},'c2_A_t',sprintf('%04g',imid),'.tif'];
        %nuclr_name{j} =['nuc/reg/xy',pos{j},'/xy',pos{j},'c3_A_t',sprintf('%04g',imid),'.tif'];
    end
    
    %store output file names for combined position images in array index n
    thrsh_name{posn+1} =['xy',outname,'c1_pha_t',sprintf('%04g',imid),'.tif'];
    phase_name{posn+1} =['xy',outname,'c1_t',sprintf('%04g',imid),'.tif'];
    fluor_name{posn+1} =['xy',outname,'c2_t',sprintf('%04g',imid),'.tif'];
    nuclr_name{posn+1} =['xy',outname,'c3_t',sprintf('%04g',imid),'.tif'];
    
    %create array for storing names of imported positional images, which will be
    %concatenated horizontally as they should be aligned already via
    %previous code
    I0 = cell(1, posn);
    Iph = cell(1, posn);
    Iflr = cell(1, posn);
    Inuc = cell(1, posn);
    
    %load the image names of all positional images
    for j=1:posn
        I0{j} = imread(thrsh_name{j});
        Iph{j} = imread(phase_name{j});
        Iflr{j} = imread(fluor_name{j});
        Inuc{j} = imread(nuclr_name{j});
    end
    
    %initialize combined image
    temp1 = [];
    temp2 = temp1;
    temp3 = temp1;
    temp4 = temp1;
    
    %for each of the positions to be concatenated; add to the right of
    %existing combined image
    for j=1:posn
        temp1 = cat(2, temp1, I0{j});
        temp2 = cat(2, temp2, Iph{j});
        temp3 = cat(2, temp3, Iflr{j});
        temp4 = cat(2, temp4, Inuc{j});
    end
    
    %export concatenated image, with xy_outname
    %exports to appropriate channel folder, subfolder stitched, and subfolder with xy_outname
    imwrite(temp1, strcat('xy',outname,'/thresh/',thrsh_name{posn+1})); 
    imwrite(temp2, strcat('xy',outname,'/phase/',phase_name{posn+1}));
    imwrite(temp3, strcat('xy',outname,'/flu/',fluor_name{posn+1}));
    imwrite(temp4, strcat('xy',outname,'/nuc/',nuclr_name{posn+1})); 
end

%clock