function flip = checkMouseSide(vid,Bkg,N)
% CHECKMOUSEISDE Checks on which side the mouse starts its run.
%
% INPUT:
% vid: MATLAB video structure or the path to a video.
% Bkg: Background image or a path to such image.
% N: Max number of images to use for the estimation.
%
% OUTPUT:
% flip: Logical value determining wether the mouse walks from left to right
% (flip == false) or right to left (flip == true).

if ischar(vid)
    vid = VideoReader(vid);
elseif isobject(vid)
    if ~strcmpi(get(vid,'Type'),'VideoReader')
        error('vid must be a VideoReader structure or the path to a video');
    end
end

if ischar(Bkg)
    Bkg = imread(Bkg);
end

if ~exist('N','var')
    N = min(10,vid.NumberOfFrames);
end

if N >= vid.NumberOfFrames
%     N = vid.NumberOfFrames;
    step = 1;
else
    step = ceil(vid.NumberOfFrames/N); % No need to be exact...
end

step_vec = 1:step:vid.NumberOfFrames;
Nsteps = length(step_vec);
C = NaN(1,Nsteps);
x_ind = 1:vid.Width;

% Checking if it is an RGB video...
sizeI = size(read(vid,1));
isrgb = length(sizeI) >2;

for i_images = 1:Nsteps
    % Computes the centroid of the white region of the image
    I = read(vid,step_vec(i_images));
    if isrgb
        I = rgb2gray(I);
    end
    
    I = im2bw(medfilt2(I - Bkg,[10 10]),0.1);
    w = sum(I,1);
    C(i_images) = sum(w.*x_ind)/sum(w);
end

% Checking if the centroid moves mostly right or left:
flip = sum(diff(C))<0;