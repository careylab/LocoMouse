function wm = weightedMean(pos_ij, C, box_size)
% WEIGHTEDMEANDETECTIONS Recomputes pos_ij as the weighted mean of the
% scores in C of a box of box_size cetered in pos_ij.

% Making sure pos_ij has the correct dimensions (2xN_candidates);
flag = false;
sizepos_ij = size(pos_ij);
[msize,msize_pos] = max(sizepos_ij);
if msize_pos == 1 || (msize_pos == 2 && sizepos_ij(1) == 1)
    pos_ij = pos_ij';
    N_candidates = 1;
    flag = true;
else
    N_candidates = msize;
end
clear msize sizepos_ij

% Zeroing negative values weights:
C(C<0) = 0;
sizeC = size(C);
% Extracting the weights:
W = getWindowFromImage(C,pos_ij,box_size);clear C
W = reshape(W,[],size(W,3));

% Computing template and local centroid:
[template_j, template_i] = meshgrid(1:box_size(1),1:box_size(2)); 
local_j = round((template_j(:)'*W)./sum(W,1)) - floor(box_size(2)/2);
local_i = round((template_i(:)'*W)./sum(W,1)) - floor(box_size(1)/2);
% transforming local coordinates into displacement vectors:
wm = pos_ij + [local_i;local_j];

if flag
    wm = wm';
end