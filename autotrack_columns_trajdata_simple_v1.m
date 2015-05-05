function autotrack_columns_trajdata_simple_v1(imN, pos, colN, flrn)



clock

outputdata =['output_xy',pos,'.mat']; 

% imN = input('get total number of images: ');
% prefix = 'xy343536'; %input('get the prefix name of image files: (string) ');
% % eg. 'xy06' or 'xy010203' or 'xy111213'
% colN = input('get the number of traps: ');
% outputdata =input('get the output data name: (string) ');
% % usually 'output_xy0x.mat' eg. 'output_xy01.mat'
% % or 'output_xy123.mat' for combined positions xy01, xy02, xy03


Traj = zeros(imN, flrn, colN);


    
fluor_name = cell(1, flrn);
Iflr = cell(1, flrn);

for imid =1:imN
    
    % get input image names
    prephase_name =['xy',pos,'/c1/xy',pos,'_c1_t',sprintf('%04g',imid),'.tif']; %c1
    for fk =1:flrn
        fluor_name{fk} =['xy',pos,'/c',num2str(fk+1),'/xy',pos,'_c',num2str(fk+1),'_t',sprintf('%04g',imid),'.tif'];
    end

    
    mask_name = ['xy',pos,'/mask/xy',pos,'ph_mask_raw_t',sprintf('%04g',imid),'.tif'];
        
    I0= imread(prephase_name);
    for fk=1:flrn
        Iflr{fk} = imread(fluor_name{fk});
    end
    
    Imask = imread(mask_name);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% initialize the column mask %%%%%%
if imid ==1

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
    
    [I3_lb colN]= bwlabeln(I3_d);
    Columns = I3_lb;
end

clear prop
prop(colN,1) = struct('c2',[],'xy',[],'center',[]);



for i=1:colN
    
    empty_tag=0; % whether the trap is empty
    
    curr_col = (Columns ==i);
    
    Ic2 = cell(1,flrn);
    for fk=1:flrn
        Ic2{fk} = double(Iflr{fk}).*curr_col;
    end
    
    
    mask = ((Imask>0).*curr_col >0);
    
    [mask_lb, nucN] = bwlabel(mask);
   
    if nucN~=0
         % fluor prop of the column
             col_c2 = cell(1, flrn);
             for fk =1:flrn
                col_c2{fk} = regionprops((mask_lb>0),Ic2{fk},'MeanIntensity','Centroid'); 
             end
             center0 = [col_c2{1}(:).Centroid];
                 
         prop(i).center = [center0(1:2:end)', center0(2:2:end)'];
         [y bott_id] = sort(prop(i).center(:,2),'descend');
         prop(i).c2 = []; 
        
        
        for c2_count =1:flrn           
            prop(i).c2 = cat(1, prop(i).c2, [col_c2{c2_count}.MeanIntensity]);
        end

        [y bott_id] = sort(prop(i).center(:,2),'descend');
        prop(i).id = bott_id(1);
        
        prop(i).btup_id = bott_id; 
        prop(i).xy = prop(i).center(bott_id(1),:);
            
    
    else
        empty_tag=1;
    end
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% assign value %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if empty_tag ==0
        curr_id = prop(i).id;
        Traj(imid, 1:flrn, i) = prop(i).c2(:,curr_id)';
    end
end

end
    
save(['uncorrected',outputdata])

clock

correct_Traj = zeros(size(Traj));
for k=1:colN
    
    k
    cutoffT = input(['get the correct start/end frames [start end] (otherwise press enter) for trap ',num2str(k),': ']);
    if ~isempty(cutoffT)
        Ts = cutoffT(1); Te = cutoffT(2);
        correct_Traj(Ts:Te,:,k) = Traj(Ts:Te,:, k);
    else
        correct_Traj(:,:,k) = Traj(:,:, k);
    end
    
    save(outputdata, 'correct_Traj');

end

end




