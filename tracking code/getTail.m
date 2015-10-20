function [Tail,T_mask] = getTail(Ii, mirror_line, box_size, w, rho, detection_thresholds,N)
% TAILDETECTION detects the tail based on thresholding.
%
% INPUT:
% I: I can be a grayscale image or a 1x2 cell. If I is a cell, it must have
% the bottom image on I{1} and the top image on I{2};
%
% mirror_line: Height at which to split the image in half if I is a
% grayscale image.
%
% box_size: A 1x2 cell with the box size for the tail on each view (1 -
% bottom, 2 - top).
%
% w: a 1x2 cell with the filter for the tail on each view (1 - bottom, 2 -
% top).
%
% rho: a 1x2 cell with the bias value for the tail filte on each view (1 -
% bottom, 2 - top).
% detection_thresholds: 
%
% N: Number of segments to divide the tail in.
%
% OUTPUT:
% Tail:
%
% T_mask: A mask to remove the tail pixels (and respective boxes) from the
% image.

switch class(Ii)
    case 'uint8'
        I_cell = cell(1,2);
        [I_cell{:}] = splitImage(Ii,mirror_line);
    case 'cell'
        I_cell = Ii;
    otherwise
        error('Supported classes for I are uint8 or cell.');
end
clear I

Tail = cell(1,2);
T_mask = cell(1,2);

if ~exist('N','var')
N = 10;
end

for i_views = 1:2
    % Detecting the tail:
    Ii = I_cell{i_views};  
%     Ii = medfilt2(Ii,[5 5]);
%     Ii = im2bw(Ii,0.1);
    T_mask{i_views} = true(size(Ii));
    
    C = conv2(double(Ii),w{i_views},'same')-rho{i_views};
    C(C < detection_thresholds(i_views)) = 0;
    BW = logical(C);
    % Finding object that maximizes the sum of scores of all pixels:
    CC = bwconncomp(BW);
%     clear BW
    fun = @(x)(sum(C(x)));
    [~,idmax] = max(cellfun(fun,CC.PixelIdxList));
    
    
    % Creating the tail mask for the current view:
    T_mask{i_views} = blackOutImageRegions(T_mask{i_views},CC.PixelIdxList{idmax}',box_size{i_views});
    
    % Estimating the tail positions:
    BW = false(size(BW));
    BW(CC.PixelIdxList{idmax}') = true;
    bb1 = sum(BW,1)>0;
    st = find(bb1,1,'first');
    en = find(bb1,1,'last');
    
    % Look for N sections measured along the x axis as it most accurately
    % represents the length of the tail:
    steps = [st st+((1:N-1).*round((en-st+1)/N)) en]; 
%     steps = [st:round((en-st+1)/N):en en];
    Ntail = length(steps)-1;
    tail = zeros(2,Ntail);
    
    
    for i_tail = 1:Ntail
        
        x = sum(BW(:,steps(i_tail):steps(i_tail+1)),1);
        y = sum(BW(:,steps(i_tail):steps(i_tail+1)),2);
        c = [x*(1:length(x))'/sum(x(:));(1:length(y))*y/sum(y(:))];
        tail(:,i_tail) = [steps(i_tail);0] + round(c);
        
    end
    Tail{i_views} = tail;
end
