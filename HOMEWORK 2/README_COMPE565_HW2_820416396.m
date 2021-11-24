% COMPE 565 Homework 2: JPEG based Image Compression
% OCTOBER 29TH 2021
% Name: Jacob Bananal
% RED ID: 820416396
% Email: jbananal2509@sdsu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M-file: COMPE565_HW2_820416396.m

% Location of output image: Flooded_house.jpg (given) 

% The COMPE_HW2_820416396.m is the code for all questions

% In order to run the code or access the .m file, you can run the 
% Matlab console or use cmd line window 

% BElOW is each of the questions implentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem 1: Encoder: (Use 4:2:0 YCbCr component image)
% Implementation 1: Image => DCT => Qunatize + Zig Zag Scan
% M-file name: COMPE565_HW2_820416396.m
% Usage: NA
% Location of output image: .\Ouptuts\
% Parameters: NA
% Other parameters here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Encoder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2: Display each band (Red, Green and Blue) of the image file (15 points)
% Implementation 1: Inverese Quantize => Inverse DCT
% M-file name: COMPE565_HW2_820416396.m
% Usage: NA
% Location of output image: .Ouptuts
% Parameters: NA
% Other parameters here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Decoder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3: ** Required by guidlines 
% • Display the Error Image (by subtracting the reconstructed image form the
% original) for the luminance image. (10 points)
% • PSNR for the luminance component of the decoded image. (10 points)
% Implementation 1: Difference Image + PSNR
% M-file name: COMPE565_HW2_820416396.m
% Usage: NA
% Location of output image: .Ouptuts
% Parameters: NA
% Other parameters here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Evaluate


SaveOutputs