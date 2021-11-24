% COMPE565 HW 2
% name: Jacob Bananal
% FALL 2021
% RED ID: 820416396
% email: jbananal2509@sdsu.edu

clear all;

% read the image 
origImg = imread('Flooded_house.jpg','jpg');
% converting to YCbCr
YCbCr = rgb2ycbcr(origImg);

figure(1);
subplot(2,2,1), subimage(origImg), title("The Original image ");
subplot(2,2,2), subimage(YCbCr), title("YCbCr");

% downsizes the chromiance by 4:2:0
y_band  = YCbCr(:,:,1);
Cb_band = YCbCr(:,:,2);
Cr_band = YCbCr(:,:,3);

totalRows = size(YCbCr,1);
totalColumns = size(YCbCr,2);

fprintf("Total rows  = %d\n",totalRows);
fprintf("Total cols  = %d\n",totalColumns);

CbCopy = Cb_band;
CrCopy = Cr_band;

% this will make all chromiance values = 0
for row = 2:2:totalRows
    for col = 2:2:totalColumns
        CbCopy(row,col) = 0;
        CrCopy(row,col) = 0;
    end
end

% removes the unwanted values from the chromiance values
Cb_420 = CbCopy(1:2:end,1:2:end);
Cr_420 = CrCopy(1:2:end,1:2:end);

figure(2)
subplot(2,2,1), subimage(Cb_band), title("The Original of the Cb band")
subplot(2,2,2), subimage(Cr_band), title("The Original of the Cr band")
subplot(2,2,3), subimage(Cb_420), title("The Cb band 4:2:0")
subplot(2,2,4), subimage(Cr_420), title("The Cr  band 4:2:0")

row_cb_420 = size(Cb_420,1);
col_cb_420 = size(Cb_420,2);

row_cr_420 = size(Cr_420,1);
col_cr_420 = size(Cr_420,2);
fprintf(" The Total rows in 4:2:0 Cb band = %d\nTotal rows in 4:2:0   cb_band = %d\n",row_cb_420,col_cb_420);
fprintf(" The Total rows in 4:2:0 Cr band = %d\nTotal rows in 4:2:0 cr_band = %d\n",row_cr_420,col_cr_420);

%takes out the 1st 8x8 block from the luminance 

y_cpy = y_band;

% block processing with DCT
dct_handle = @dct2;
% 8x8 block processing on y_cpy
y_dct = blkproc(y_cpy,[8 8],dct_handle);
% now computing the DC the coefficients for the cb and cr components
cb_cpy = Cb_420;
cr_cpy = Cr_420;

cb_dct = blkproc(cb_cpy,[8 8],dct_handle);
cr_dct = blkproc(cr_cpy,[8 8],dct_handle);

figure(3);
subplot(2,2,1), subimage(y_dct), title("dct y band")
subplot(2,2,2), subimage(cb_dct), title("dct cb band")
subplot(2,2,3), subimage(cr_dct), title("dct cr band")

% taking out the 1st 2 block from 6th row dct
y_dct_block1 = y_dct(41:48,1:8);
y_dct_block2 = y_dct(41:48,9:16);

% displays the DCT coefficients of luminance blocks
fprintf("\nThe DCT coeff of 1st block:\n")
disp(y_dct_block1)
fprintf("\nThe DCT coeff of 2nd block:\n")
disp(y_dct_block2)

figure(4);
subplot(2,2,1), subimage(y_dct_block1), title("8x8 y band block1")
subplot(2,2,2), subimage(y_dct_block2), title("8x8 y band block2")

% trucating the 8x8 blocks for correct representation using imshow()
% truncate function defined at the end before the zig zag 
t_blk1 = truncate(y_dct_block1);
fprintf("The Truncated 8x8 blk1:\n");
disp(t_blk1);
t_blk2 = truncate(y_dct_block2);
fprintf("The Truncated 8x8 blk2:\n");
disp(t_blk2);

% displays the truncated blocks
figure(5);
subplot(2,2,1), imshow(t_blk1),title("The truncated 8x8 block1")
subplot(2,2,2), imshow(y_dct_block1),title("The 0riginal 8x8 block1")
subplot(2,2,3), imshow(t_blk2),title("The truncated 8x8 block2")
subplot(2,2,4), imshow(y_dct_block2),title("org 8x8 block2")


% Q2: Quantization
%  quantizing the luminance and chromiance matrix with the quantizer
% matrix from the lecture.


% quantization matrix for luminance from lecture powerpoint
luminance_matrix = [16 11 10 16 24 40 51 61;12 12 14 19 26 58 60 55;14 13 16 24 40 57 69 56;
14 17 22 29 51 87 89 62;18 22 37 56 68 109 103 77;24 35 55 64 81 104 113 92;
49 64 78 87 108 121 120 101;72 92 95 98 112 100 103 99];

% quantization marix for chromiance from lecture powerpoint
chromiance_matrix = [17 18 24 47 99 99 99 99;18 21 26 66 99 99 99 99;24 26 56 99 99 99 99 99;
47 66 99 99 99 99 99 99;99 99 99 99 99 99 99 99;99 99 99 99 99 99 99 99;
99 99 99 99 99 99 99 99;99 99 99 99 99 99 99 99];

% create copies of y,cb and cr dct matricies
Y_DCT_cpy  = y_dct;
Cb_DCT_cpy = cb_dct;
Cr_DCT_cpy = cr_dct;

% divides the dct block structure with the lum matrix.
q_lum = @(Y_DCT_cpy)round(Y_DCT_cpy./luminance_matrix);
% 8x8 block processing on the quantized matrix
Y_DCT_q = blkproc(Y_DCT_cpy,[8 8],q_lum);

% dividing the dct block_struct data with the chromiance matrix.
q_chr1 = @(Cb_DCT_cpy)round(Cb_DCT_cpy./chromiance_matrix);
q_chr2 = @(Cr_DCT_cpy)round(Cr_DCT_cpy./chromiance_matrix);
% 8x8 block processing on the quantized matrix
Cb_q_DCT = blkproc(Cb_DCT_cpy,[8 8],q_chr1);
Cr_q_DCT = blkproc(Cr_DCT_cpy,[8 8],q_chr2);

% reports the data for the first 2 blocks in the 6th row from the top.
q_y_blk1  = Y_DCT_q(41:48,1:8);
q_y_blk2  = Y_DCT_q(41:48,9:16);
q_cb_blk1 = Cb_q_DCT(41:48,1:8);
q_cb_blk2 = Cb_q_DCT(41:48,9:16);
q_cr_blk1 = Cr_q_DCT(41:48,1:8);
q_cr_blk2 = Cr_q_DCT(41:48,9:16);

% displays the quantized y,cb,cr blocks
figure(6);
subplot(3,2,1), subimage(q_y_blk1), title("The quantized 1st 8x8 block of y band")
subplot(3,2,2), subimage(q_y_blk2), title(" The quantized 2nd 8x8 block of y band")
subplot(3,2,3), subimage(q_cb_blk1), title("The quantized 1st 8x8 block of cb band")
subplot(3,2,4), subimage(q_cb_blk2), title("The quantized 2nd 8x8 block of cb band")
subplot(3,2,5), subimage(q_cr_blk1), title("The quantized 1st 8x8 block of cr band")
subplot(3,2,6), subimage(q_cr_blk2), title("The quantized 2nd 8x8 block of cr band")

% prints out the DC coeff of the luminance block
fprintf("The DC coeff for 1st 8x8 block of Y band:%d\nThe DC coeff for 2n 8x8 block of Y band:%d\n",q_y_blk1(1,1),q_y_blk2(1,1));
fprintf("The 8x8 luminance block1:\n");
disp(q_y_blk1);


% finds the AC coeffs
% uses the zizag scanning to find the AC cefficients.
% zigzag() function defined at the end of this script.
ac_mat1 = zigzag(q_y_blk1);
fprintf("The Ac coeffs of luminance component 8x8 blk1:\n");
display(ac_mat1(2:end));
ac_mat2 = zigzag(q_y_blk2);
fprintf("The 8x8 luminance block2:\n");
disp(q_y_blk2);
fprintf("\nThe Ac coeffs of luminance component 8x8 blk2:\n");
display(ac_mat2(2:end));



% Decoder Part:

% Q3: Inverse Quantization
% inverse quantize the quantized image


% mutiplies the quantized lum matrix with the quantization matrix for
% inverse dct operation

y_q_cpy  = Y_DCT_q;
cb_q_cpy = Cb_q_DCT;
cr_q_cpy = Cr_q_DCT;

% implements the inverse quantization by multiplying the quantization
% matrix to get the inverse quantized matrices
iq_lum = @(y_q_cpy)round(y_q_cpy.*luminance_matrix);
% 8x8 block processing of the entire matrix
iq_y = blkproc(Y_DCT_q,[8 8],iq_lum);

% performs the same operation for cb and cr
iq_chr1 = @(cb_q_cpy)round(cb_q_cpy.*chromiance_matrix);
iq_chr2 = @(cr_q_cpy)round(cr_q_cpy.*chromiance_matrix);

iq_cb = blkproc(Cb_q_DCT,[8 8],iq_chr1);
iq_cr = blkproc(Cr_q_DCT,[8 8],iq_chr2);

% Displays the inverse quantized y,cb and cr images.
figure(7);
subplot(3,1,1), subimage(iq_y), title("The inverse quantized image of y band")
subplot(3,1,2), subimage(iq_cb), title("The inverse quantized image of cb band")
subplot(3,1,3), subimage(iq_cr), title("The inverse quantized image of cr band")


% Q4 reconstruct, psr, error image of y band


% now will apply the inverse DCT to y cb and cr component matrices
% separately.

idct_handle = @(block_struct)idct2(block_struct.data);
% using the above handle with the blockproc()
y_idct  = blockproc(iq_y,[8 8],idct_handle);
cb_idct = blockproc(iq_cb,[8 8],idct_handle);
cr_idct = blockproc(iq_cr,[8 8],idct_handle);

% since ycbcr needs to be in the unsigned int format typecasting the above
% array to uint16

y_idct  = uint8(y_idct);
cb_idct = uint8(cb_idct);
cr_idct = uint8(cr_idct);

% displays the inverse DCT y,cb and cr bands
figure(8)
subplot(3,1,1), subimage(y_idct), title("The inverse DCT image of y band")
subplot(3,1,2), subimage(cb_idct), title("The inverse DCT image of cb band")
subplot(3,1,3), subimage(cr_idct), title("The inverse DCT image of cr band")


% using linear interpolation method
t_r = size(YCbCr,1);
t_c = size(YCbCr,2);
rep_ycbcr = zeros(t_r,t_c,3);
y_idct_cpy = y_idct;
cb_idct_cpy = cb_idct;
cr_idct_cpy = cr_idct;

rep_ycbcr(1:2:t_r,1:2:t_c,2) = cb_idct_cpy(:,:);
rep_ycbcr(1:2:t_r,1:2:t_c,3) = cr_idct_cpy(:,:);
% temp var to hold upsampled cb and cr values
temp_cb = rep_ycbcr(:,:,2)
temp_cr = rep_ycbcr(:,:,3);


% performs linear interpolation on cb and cr bands for upsampling
if mod(t_r,2) == 0
    for r = 2:2:t_r - 2
        for c = 1:2:t_c
            temp_cb(r,c) = round((temp_cb(r-1,c)+temp_cb(r+1,c))/2);
            temp_cr(r,c) = round((temp_cr(r-1,c)+temp_cr(r+1,c))/2);
        end
    end
    for c = 2:2:t_c - 2
        for r = 1:1:t_r
            temp_cb(r,c) = round((temp_cb(r,c-1)+temp_cb(r,c+1))/2);
            temp_cr(r,c) = round((temp_cr(r,c-1)+temp_cr(r,c+1))/2);
        end
    end
    for r = t_r
        for c = t_c
            temp_cb(r,c) = temp_cb(r-1,c);
            temp_cr(r,c) = temp_cr(r-1,c);
        end
    end
else
    for r = 2:2:t_r - 1
        for c = 1:2:t_c
            temp_cb(r,c) = round((temp_cb(r-1,c)+temp_cb(r+1,c))/2);
            temp_cr(r,c) = round((temp_cr(r-1,c)+temp_cr(r+1,c))/2);
        end
    end
    for c = 2:2:t_c - 1
        for r = 1:2:t_r
            temp_cb(r,c) = round((temp_cb(r,c-1)+temp_cb(r,c+1))/2);
            temp_cr(r,c) = round((temp_cr(r,c-1)+temp_cr(r,c+1))/2);
        end
    end
end

rep_ycbcr = cat(3,y_idct,temp_cb,temp_cr);
% converts the reconstructed YCbCr image to RGB
rep_rgb = ycbcr2rgb(rep_ycbcr);

figure(9);
subplot(2,2,1), subimage(YCbCr), title("Original YCbCr Image")
subplot(2,2,2), subimage(rep_ycbcr), title("The Reconstructed YCbCr image")
subplot(2,2,3), subimage(origImg), title("Original RGB image")
subplot(2,2,4), subimage(rep_rgb), title("The Reconstructed RGB image")


% Calculates the error in the image only for lumiance.

% subtracing the reconstructed y band from original y band
error = YCbCr(:,:,1) - rep_ycbcr(:,:,1);

% displays the error image
figure(10);
subplot(3,1,1), subimage(YCbCr(:,:,1)), title("The original y band")
subplot(3,1,2), subimage(rep_ycbcr(:,:,1)), title("The reconstructed y band")
subplot(3,1,3), subimage(error), title("The error between reconstructed and original image")


% using the psnr() matlab function to calucate the psnr value.
psnr_lum = psnr(rep_ycbcr(:,:,1),YCbCr(:,:,1));
fprintf("The Peak SNR of decoded Y band:%f\n",psnr_lum);

%The Peak SNR of decoded Y band:33.519912


% this defines a function to truncate 8x8 blocks for imshow()

function op = truncate(blk)
min_v = min(min(blk));
min_v = -min_v;
op = blk;
for i = 1:1:8
    for j = 1:1:8
        if op(i,j) == 0
            op(i,j) = op(i,j);
        else
            op(i,j) = op(i,j) + min_v;
        end
    end
end
op = op./255;
end



% zigzag function 

function AC_matrices = zigzag(q_blk)
index_point=1;
mat_dct = q_blk;
AC_matrices = [];
N = 8;
for c=1:2*N-1
    if c<=N
        if mod(c,2)==0
            j=c;
            for i=1:c
                AC_matrices(index_point)= mat_dct(i,j);
                index_point=index_point+1;j=j-1;
            end
        else
            i=c;
            for j=1:c
                AC_matrices(index_point)= mat_dct(i,j);
                index_point=index_point+1;i=i-1;
            end
        end
    else
        if mod(c,2)==0
            p=mod(c,N); j=N;
            for i=p+1:N
                AC_matrices(index_point)=mat_dct(i,j);
                index_point=index_point+1;j=j-1;
            end
        else
            p=mod(c,N);i=N;
            for j=p+1:N
                AC_matrices(index_point)=mat_dct(i,j);
                index_point=index_point+1;i=i-1;
            end
        end
    end
end
end