% [nouse clump_idx0] = find(PAratio > 1.5);
% clump_idx= setdiff(clump_idx0, edgeid);


function [PAratio_out, mask_out, mask_label_out] = de_clump_v2(mask_in, label_mask, PAratio, clump_idx)

I_cell_lb = label_mask;
I_cell = mask_in;

% [nouse pre_clump_idx] = find((PAratio > 1.35).*(PAratio <1.5)); 
% 
%     % for these <1.4<1.5, test bwconvhull
%     if ~isempty(pre_clump_idx)
%         for kk = 1:length(pre_clump_idx)
%             clump_mask0 = (I_cell_lb == pre_clump_idx(kk));
%             clump_mask1 = bwconvhull(clump_mask0);
%             clump_prop0 = regionprops(clump_mask1, 'Area','Perimeter');
%             
%             clump_ratio0 = [clump_prop0(:).Perimeter].^2./[clump_prop0(:).Area]/4/pi;
%             if clump_ratio0 <1.35
%                     
%                     I_cell = I_cell -clump_mask0 + clump_mask1;
%             end
%         end
%     end
            
        
   
    
    
    if ~isempty(clump_idx)
        for kk = 1:length(clump_idx)
            clump_mask0 = (I_cell_lb == clump_idx(kk));
            
            label_ini =1; 
            label_curr = 1;
            erode_count =0;
            fixed_tag =0;
            
            PAratio_clump = PAratio(clump_idx(kk));
                        
                        
            % test imdilate disk 1 first see if the ratio drop
            % in case some
            clump_ini = imdilate(clump_mask0, strel('disk',1));
            clump_prop_ini = regionprops(clump_ini, 'Area','Perimeter');
            
            clump_ratio_ini = [clump_prop_ini(:).Perimeter].^2./[clump_prop_ini(:).Area]/4/pi;
            
            if clump_ratio_ini <1.5
                PAratio_clump = clump_ratio_ini;
                clump_mask2 = bwconvhull(clump_mask0);
                I_cell = I_cell - clump_mask0 + clump_mask2;

            else        
            % use dilation to separate cell clumps.
%                     clump_ini = imdilate(clump_mask0, strel('disk',2));
                    clump_mask = bwmorph(clump_mask0,'thin',1);
                    clump_prop = regionprops(clump_mask, 'Area','Perimeter');
                    label_curr = size(clump_prop,1);
                    clump_ratio = [clump_prop(:).Perimeter].^2./[clump_prop(:).Area]/4/pi;
                    PAratio_clump = max(clump_ratio);
                    
                    while label_curr<2 && PAratio_clump > 1.5
                    % label_curr =1, PAratio_clump >1.5  keep
                    % label_curr =1, PAratio_clump <1.5  jump
                    % label_curr >1, PAratio_clump any   jump
                        clump_mask2 = imerode(clump_mask, strel('disk',1));

                        clump_prop = regionprops(clump_mask2, 'Area','Perimeter');

                        label_curr = size(clump_prop,1);
                        clump_ratio = [clump_prop(:).Perimeter].^2./[clump_prop(:).Area]/4/pi;
                        PAratio_clump = max(clump_ratio);
                        
%                         clump_tiny = find([clump_prop.Area]<5);
%                         
%                         if ~isempty(clump_tiny)
%                            label_curr=length([clump_prop.Area]) - length(clump_tiny);
%                            clump_big = find([clump_prop.Area]>4);
%                            clump_big_lb = bwlabeln(clump_mask2);
%                            
%                            for jj =1:length(clump_big)
%                                 clump_mask2 = (clump_big_lb==clump_big(1))+(clump_big_lb==clump_big(jj));
%                            end
%                            clump_mask2 = clump_mask2>0;  erode_count = 1;
%                         end
                        
                        if isempty(clump_ratio)
                            break;
                        end

                        PAratio_clump = max(clump_ratio);
                        clump_mask = clump_mask2;
                        erode_count = erode_count+1;

                    end
                    
                    if erode_count>4
                        
%                         figure; imagesc(clump_ini)
                        clump_ini = imdilate(clump_mask0, strel('disk',2));
                        clump_prop_ini = regionprops(clump_ini, 'Area','Perimeter');

                        clump_ratio_ini = [clump_prop_ini(:).Perimeter].^2./[clump_prop_ini(:).Area]/4/pi;

                        if clump_ratio_ini <1.5
                            PAratio_clump = clump_ratio_ini;
                            clump_mask2 = bwconvhull(clump_mask0);
                            I_cell = I_cell - clump_mask0 + clump_mask2;
                            fixed_tag=1;
                        end
                        
                        
                    end
                        
                        
                        
                    % finish while loop
                    
                    % if it's divided into two objects
                    if fixed_tag==0
                        if label_curr>1 

                            clump_lb = bwlabel(clump_mask);

                            expand_mask = zeros([size(clump_lb),label_curr]);
                            for jj =1:label_curr
                                expand_mask(:,:,jj) = imdilate((clump_lb==jj),strel('disk',erode_count));
                            end
                            % if connected, use line not disk
                            check_conn= bwconncomp(sum(expand_mask,3),4);
                            temp = (sum(expand_mask,3)>0);
                            
                            if check_conn.NumObjects==1 
                                % didn't separate them
                                
                                temp2 = imerode(temp, strel('disk',1));
                                I_cell = (I_cell - clump_mask0 + temp2)>0;
                            else
                                %two separated objs
                                % get the "neck" region and dilate
                                I_cell = (I_cell - clump_mask0 + temp)>0;
%                                 remove_neck0 = clump_mask0 - sum(expand_mask,3);
%                                 remove_neck = imdilate(remove_neck0, strel('line',3,0));
%                                 I_cell = (I_cell- remove_neck)>0;
                            end
                            
                            

                        elseif isempty(PAratio_clump)
                            I_cell = I_cell - clump_mask0;

                        elseif PAratio_clump<1.5
                            expand_mask = bwconvhull(clump_mask0);
                            % instead of
                            % expand_mask = imdilate(clump_mask,strel('disk',erode_count));
                            I_cell = I_cell- clump_mask0 + expand_mask;
                        % if label_curr=1 and PAratio_clump <1.5? => one single  cell
                        end
                        
                        fixed_tag=1;
                    end
            end
            
        end
    end
    
    mask_out = (I_cell)>0;
    mask_label_out = bwlabel(mask_out);
    
    prop_out = regionprops(mask_out, 'Area','Perimeter');

    PAratio_out = [prop_out(:).Perimeter].^2./[prop_out(:).Area]/4/pi;
 
    
end
    
    