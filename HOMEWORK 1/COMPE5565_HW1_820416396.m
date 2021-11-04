% Compe565 homework 1 
% name: Jacob Bananal
% RED ID: 820416396
% email: jbananal2509@sdsu.edu

clear all;

% problem 1: Read & display the image 

origImg = imread('Flooded_house.jpg','jpg');
figure(1);
imshow(origImg),title('Original Image'), axis on;
totalRows = size(origImg,1);
totalColumns = size(origImg,2);

% problem 2: Displays the RBG bands of the given original image 

% makes 3 copies of the original image
RedBand = origImg;
BlueBand = origImg;
GreenBand = origImg;

% red component
RedBand(:,:,2)=0;
RedBand(:,:,3)=0;
% green component 
GreenBand(:,:,1)=0;
GreenBand(:,:,3)=0;
% blue component
BlueBand(:,:,1)=0;
BlueBand(:,:,2)=0;

figure(2);
subplot(2,2,1), imshow(origImg), title('Original Image'), axis on;
subplot(2,2,2), imshow(RedBand), title('Red Band'), axis on;
subplot(2,2,3), imshow(GreenBand), title('Green Band'), axis on;
subplot(2,2,4), imshow(BlueBand), title('Blue Band'), axis on;

% problem 3: 

% part 3.1: convert an RGB image into a YCbCr image
YCbCr = rgb2ycbcr(origImg); 
% part 3.2: convert a YCbCr image into RGB format
RGB = ycbcr2rgb(YCbCr); 

figure(3);
subplot(1,2,1),imshow(YCbCr),title('RGB to YCbCr');
subplot(1,2,2),imshow(RGB),title('YCbCr toRGB');

% problem 4: Displays each Y band, Cb band, and Cr band seperately 

y_band = YCbCr(:,:,1);
cb_band = YCbCr(:,:,2);
cr_band = YCbCr(:,:,3);
figure(4);
subplot(2,2,1), imshow(YCbCr), title('YCbCr'), axis on;
subplot(2,2,2), imshow(y_band), title('Y Band'), axis on;
subplot(2,2,3), imshow(cb_band), title('Cb Band'), axis on;
subplot(2,2,4), imshow(cr_band), title('Cr Band'), axis on;

% problem 5: subsample Cb and Cr bands using 4:2:0 and display both bands
CbSub = cb_band;
CrSub = cr_band;
 % goes through a FOR-statement
for row = 2:2:totalRows
    for col = 2:2:totalColumns
        CbSub(row, col)=0;
        CrSub(row, col)=0;
    end
end

Cb_420 = CbSub(1:2:end, 1:2:end);
Cr_420 = CrSub(1:2:end, 1:2:end); 
figure(5);
subplot(1,2,1), imshow(Cb_420), title('4:2:0 Subsample of Cb band'), axis on;
subplot(1,2,2), imshow(Cr_420), title('4:2:0 Subsample of Cr band'), axis on;

% problem 6: Upsample and display the Cb and Cr bands using:
% 6.1 Linear interpolation
li_cb = cb_band;
li_cr = cr_band;
     
for row = 1:totalRows-1
    for col = totalColumns-1
        if (li_cb(row, col) == 0 && li_cr(row, col) == 0)
            if(mod(row,2)==0)
                li_cb(row, col) = (li_cb(row-1, col)/2) + (li_cb(row+1, col)/2); 
                li_cr(row, col) = (li_cr(row-1, col)/2) + (li_cr(row+1, col)/2); 
            end
        end   
    end
end

for row = 1:totalRows-1
    for col = 1:totalColumns-1
        if (li_cb(row, col) == 0 && li_cr(row, col) == 0) 
            if(mod(col,2)==0)   
                li_cb(row, col) = (li_cb(row, col-1)/2) + (li_cb(row, col+1)/2); 
                li_cr(row, col) = (li_cr(row, col-1)/2) + (li_cr(row, col+1)/2); 
            end
        end
        
    end
end

figure(6);
subplot(1,2,1); imshow(li_cb); title('Linear interpolation of Cb band'), axis on;
subplot(1,2,2); imshow(li_cr); title('Linear interpolation of Cr band'), axis on;

%6.2 Simple row or column replication.
rcr_cb = cb_band;
rcr_cr = cr_band;
for row = 1:totalRows
    for col = 1:totalColumns
        if (rcr_cb(row, col) == 0 && rcr_cr(row, col) == 0) 
            if(mod(row,2) == 1) 
                rcr_cb(row, col) = rcr_cb(row, col-1);
                rcr_cr(row, col) = rcr_cr(row, col-1);
            end
        end
    end
end

for row = 1:536
    for col = 1:704
        if (rcr_cb(row, col) == 0 && rcr_cr(row, col) == 0)
            if(mod(row,2) == 0) 
                rcr_cb(row, col) = rcr_cb(row-1, col);
                rcr_cr(row, col) = rcr_cr(row-1, col);
            end
        end
    end
end
figure(7);
subplot(1,2,1); imshow(rcr_cb); title('Row&Column Replication of Cb band'), axis on;
subplot(1,2,2); imshow(rcr_cr); title('Row&Column Replication of Cr band'), axis on;

% problem 7: Covert image into RGB format 

LI_total = cat(3, y_band, li_cb, li_cr);
RC_total = cat(3, y_band, rcr_cb, rcr_cr);
LI2RGB = ycbcr2rgb(LI_total);
RC2RGB = ycbcr2rgb(RC_total);

% problem 8: Displays the original and reconstructed images 

figure(8);
subplot(1,3,1), imshow(origImg), title('Original RGB Image'), axis on;
subplot(1,3,2), imshow(LI2RGB), title('Linear Interpolation yCbCr to RGB'), axis on;
subplot(1,3,3), imshow(RC2RGB), title('RC Replication yCbCr to RGB'), axis on;

% problem 9 
% In the reconstructed file, from my observation
% while using linear interpolation, the image unsampled has a lot better
% quality than the image unsample from the Row-Column replication. This 
% happens because the row column replication pixels can be very larged
% due to the n-1 extending to n rows or columns. That is why linear 
% interpolation can be better since it only substitutes the neighboring 
% column or row into the missing pixels of the image. 

% problem 10: Measurements of the mean squared error between the original
% and the reconstructed image

pixdiff = origImg - LI2RGB;
mse = (sum(sum(pixdiff.^2)))/(totalRows*totalColumns);
fprintf('MSE of Red Band: %f \n', mse(:,:,1));
fprintf('MSE of Green Band: %f \n', mse(:,:,2));
fprintf('MSE of Blue Band: %f \n', mse(:,:,3));

%MSE of Red Band: 0.166153 
%MSE of Green Band: 0.121202 
%MSE of Blue Band: 0.197822
% there seems to be a slight distortion in the green band whereas there
% seems to be more of a distortion in the blue band 

%  problem 11: Compression Ration obtained by subsampling Cr and Cb 
% components for the 4:2:0 approach 

original_house_size = size(origImg,1)*size(origImg,2)*3;
subsample_house_size = size(y_band,1)*size(y_band,2)+size(Cb_420,1)*size(Cb_420,2)+size(Cr_420,1)+size(Cr_420,2);

ratio = original_house_size/subsample_house_size;

fprintf('Compression Ratio: %f \n', ratio);
fprintf('\nThe File size of the original YCbCr img: %d',(original_house_size)/8);
fprintf('\nThe File size of the Subsampled img: %d',(subsample_house_size)/8);

% The Compression Ratio: 2.396849 
% the compression is generated a 2:1 ratio since ther file size of the
% of the orignal and subsampled are divided to each other. 