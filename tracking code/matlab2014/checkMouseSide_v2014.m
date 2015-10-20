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

%%% Conflicts with some MATLAB versions that changed the behaviour of Video
%%% Reader.
% if ischar(vid)
%     vid = VideoReader(vid);
% elseif isobject(vid)
%     if ~strcmpi(get(vid,'Type'),'VideoReader')
%         error('vid must be a VideoReader structure or the path to a video');
%     end
% end

if ischar(Bkg)
    Bkg = imread(Bkg);
end

N_frames = vid.Duration * vid.frameRate;

if ~exist('N','var')
    N = min(10,N_frames);
end

if N >= N_frames
%     N = N_frames;
    step = 1;
else
    step = ceil(N_frames/N); % No need to be exact...
end

step_vec = 1:step:N_frames;
Nsteps = length(step_vec);
C = NaN(1,Nsteps);
x_ind = 1:vid.Width;

% Checking if it is an RGB video...
sizeI = size(readFrame(vid));
% sizeI = size(read(vid,1));
isrgb = length(sizeI) >2;

for i_images = 1:Nsteps
    % Computes the centroid of the white region of the image
    set(vid,'CurrentTime',(step_vec(i_images)-1)/vid.frameRate);
    I = readFrame(vid);
    
%     I = read(vid,step_vec(i_images));
    if isrgb
        I = rgb2gray(I);
    end
    
    I = im2bw(medfilt2(I - Bkg,[10 10]),0.1);
    w = sum(I,1);
    C(i_images) = sum(w.*x_ind)/sum(w);
end

% Checking if the centroid moves mostly right or left:
flip = sum(diff(C))<0;