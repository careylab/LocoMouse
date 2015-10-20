function I = imreadGrayscaleBrightnessAdjustment(im_path)
% IMREADGRAYSCALEBRIGHTNESSADJUSTMENT Reads an image from a file path and
% adjusts the brightness so it fills the whole range of the image (0 to 
% 255).
%
% Adjustment is done using the sc function.

I = sc(imread(im_path),'gray');I = rgb2gray(uint8(I*255));