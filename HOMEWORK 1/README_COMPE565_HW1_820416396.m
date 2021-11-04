% COMPE 565 Homework 1: Basic Digital Image Processing
% Sept 22. 2021
% Name: Jacob Bananal
% RED ID: 820416396
% Email: jbananal2509@sdsu.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% M-file: COMPE565_HW1_820416396.m

% Location of output image: Flooded_house.jpg (given) 

% The COMPE_HW1_820416396.m is the code for all 11 questions

% Question 9, Question 10, & Question 11 answers and comments 
% are in the comments section of the .m file code 

% In order to run the code or access the .m file, you can run the 
% Matlab console or use cmd line window 

% BElOW is each of the questions (1-11) implentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% QUESTIONS
% 1. Read and display the image using Matlab 
% 2. Display each band (Red, Green and Blue) of the image file 
    %rgbt

% 3. Convert the image into YCbCr color space
% 4. Display each band separately (Y, Cb and Cr bands)
    %ycbcrt


% 5. Subsample Cb and Cr bands using 4:2:0 and display both bands. 
% 6. Upsample and display the Cb and Cr bands using
% 6.1. Linear interpolation
% 6.2. Simple row or column replication
% 7. Convert the image into RGB format
% 8. Display the original and reconstructed images (the image restored from
% the YCbCr coordinate)
% 9. Comment on the visual quality of the reconstructed image for both the 
% upsampling cases
% 10. Measure MSE between the original and reconstructed images (obtained 
% using linear interpolation only)
% comments in the comments section of the code 
% 11. Comment on the compression ratio achieved by subsampling Cb and Cr 
% components for 4:2:0 approach. 
% comments in the comments section of the code
