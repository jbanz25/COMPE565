% COMPE565 Homework 3 - motion estimation for video compression
% Author: Jacob Bananal
% RED ID: 820416396
% email : jbananal2509@sdsu.edu
% FALL 2021 

clear all;
clc;

% Reading the walk video sequence.
video_obj = VideoReader('walk_qcif.avi');
i = 0;
j = 1;
v_width = video_obj.Width;
v_height = video_obj.Height;
fprintf("The dimensions of the given Video sequence:[Height x Width] :[%d x %d]\n",v_height,v_width);

% this extracts every frame from the video sequence that was given
% then converts each extracted frame to ycbcr.
while hasFrame(video_obj)
    v_frm(j).cdata = readFrame(video_obj);
    ycbcr = rgb2ycbcr(v_frm(j).cdata(:,:,:));
    y(:,:,j) = ycbcr(:,:,1);
    cb(:,:,j) = ycbcr(1:2:end,1:2:end,2);
    cr(:,:,j) = ycbcr(1:2:end,1:2:end,3);
    % this calcs the # of frames in the video sequence
    i = i+1;
    j = j+1;
end
% prints the total frames in the given video sequence 
fprintf("The total Frames in the video sequence:%d\n",i);

% this defines the Frame 6  which is the first frame from given GoP(6:10)[IPPPP] as
% intra-frame
rf = y(:,:,6);
figure(1)
imshow(rf)
title('Intra-Frame')

% this computes the total Macroblocks for the video sequence
total_macroblocks = (ceil(size(rf,1)/16))*(ceil(size(rf,2)/16));
fprintf("The total Macroblocks of the video sequence:%d\n",total_macroblocks);

% Variance which will compute the total comparisons done using algorithm
% below
compar = 0;
% this will go through a the loop to find difference matrix ,motion estimation
% and motion vector, therefore this will start at the require 6 then to 10
for i = 6:9
    % this extracts the Y,cb and Cr compoments, then downsamples the chromiance
    % of the video frame img.
    y(:,:,i-1)  = ycbcr(:,:,1);
    cb(:,:,i-1) = ycbcr(1:2:end,1:2:end,2);
    cr(:,:,i-1) = ycbcr(1:2:end,1:2:end,3);
    
    % setting the current and reference frames for motion prediction
    target_frm = y(:,:,i);
    refence_frm = y(:,:,i+1);
    target_frm_2 = double(target_frm);
    reference_frm_2 = double(refence_frm);
    
    % defining the MB size as 16 for Y comp.
    macroblock_size = 16;
    
    % getting total rows and column in the current frame
    [row,col] = size(target_frm_2);
    temp_diff_frm = zeros(macroblock_size,macroblock_size);
    disp_vect = zeros(1,2);
    
    % two matrices for holding the diff frame and tgt frame
    search_window = zeros(total_macroblocks,2,2);
    % var to compute the search window movement.
    c = 1;
    for r_macroblock = 1:macroblock_size:row
        for c_macroblock = 1:macroblock_size:col
            % this extracts the temporary target frame window.
            temp_target_macroblock = target_frm_2([r_macroblock:r_macroblock+macroblock_size-1],[c_macroblock:c_macroblock+macroblock_size-1]);
            max_val = 66000;
            for blk_r = -8:8
                for blk_c = -8:8
                    % this will pad the extra rows and columns inside the
                    % search window
                    inc_row = r_macroblock+blk_r;
                    inc_col = c_macroblock+blk_c;
                    if ((inc_row + 16 - 1) <= row) && ((inc_col +16 - 1) <= col) && (inc_row > 0)  && (inc_col > 0)
                        ref_macroblock_sw = reference_frm_2(inc_row:inc_row+macroblock_size-1,inc_col:inc_col+macroblock_size-1);
                        temp_diff_frm = temp_target_macroblock - ref_macroblock_sw;
                        
                        % computes the MSE for current Macroblock and
                        % reference Macroblock
                        % from the computed search window
                        mse_macroblock = sum(sum(temp_diff_frm.^2));
                        mse_macroblock = mse_macroblock./256;
                        % this finds the closest block matching which will
                        % also compute
                        % the MAD
                        if mse_macroblock < max_val
                            max_val = mse_macroblock;
                            disp_vect = [inc_row - r_macroblock,inc_col - c_macroblock];
                            rc_img(r_macroblock:r_macroblock+macroblock_size -1,c_macroblock:c_macroblock+macroblock_size -1) = reference_frm_2(inc_row:inc_row+macroblock_size -1,inc_col:inc_col+macroblock_size -1);
                          
                            % increments the compare count to display
                            % total MAD comaprisons at the end of the loop
                            compar = compar +1;
                        elseif mse_macroblock == max_val
                            pad_r_c = (r_macroblock - inc_row)^2 + (c_macroblock - inc_col)^2;
                            if pad_r_c < max_val
                                disp_vect = [inc_row - r_macroblock,inc_col - c_macroblock];
                            end
                        end
                    end
                end
            end
            diff_frm(r_macroblock:r_macroblock+macroblock_size -1,c_macroblock:c_macroblock+macroblock_size -1) = temp_diff_frm;
           
            % this will then update the  the search window
            search_window(c,:,1) = [r_macroblock,c_macroblock];
            search_window(c,:,2) = disp_vect;
            c = c+1;
        end
    end
    fprintf("\nThis is the MSE of Frame:%d and Frame:%d\n",(i),(i+1));
    disp(mse_macroblock);
    figure()
    % displays the Motion Vector Representation using the quiver command as
    % mentioned in the required 
    quiver(search_window(:,2,1), search_window(:,1,1), search_window(:,2,2), search_window(:,1,2));
    title(['The Motion vector Representation of Y component image frame:',num2str(i),'and image frame:',num2str(i+1)]);
    grid on
    reconst_img = uint8(rc_img);
    figure()
    subplot(2,2,1),imshow(target_frm),title(['Y component Original Video Frame:',num2str(i)])
    subplot(2,2,2),imshow(reconst_img),title(['Predicted Target Y component of the Video Frame:',num2str(i+1)])
    pred_err = target_frm - reconst_img;
    subplot(2,2,3),imshow(pred_err),title('The error between the predicted and target frame')
end
[add,sub] = exhaustive_search_load_cal(16);
fprintf("The total additions for while computing:%d\nTotal substractions for while computing:%d\nTotal comparisons for while computing:%d\n",add,sub,compar);


% exhaustive_search_load_cal()
% returns the total additions and substractions performed
% during computation of the motion estimation which will be called at the
% end of the code above. 
function [total_add,total_sub] = exhaustive_search_load_cal(macroblock_size)
    total_s = ((2*8)+1)^2;
    total_add = (2*(macroblock_size^2))*total_s;
    total_sub = (macroblock_size^2)*total_s;
end
