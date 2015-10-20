function Ibw = blackOutImageRegions(Ibw, centers, box_size)
% BLACKOUTIMAGEREGIONS  masks out image regions given locations and box
% size.
%
% Input:
%
% Ibw: an IxJ logical matrix.
% centers: 1xN or 2xN indices reprensenting image locations.
% box_size: the region around the centers to be masked out (set to false).
% canonical_row_col: The canonical row and column indices that must be
% added to the provided center coordinates to result on the box on the
% appropriate region, provided as a 1x2 cell where canonical_row_col{1}
% corresponds to the canonical row and canonical_row_col{2} to the column.
%
% Output:
% 
% Ibw: The result from masking regions form the input image.
% canonical_row_col: If the canonical_row_col variable has not been
% provided, it is computed and returned to be used in future
% computations.

size_I = size(Ibw);
Nc1 = size(centers,1);
if Nc1 == 1
    % Convert 2D coordinates to linear indices:
    c = cell(2,1);
    [c{:}] = ind2sub(size_I,centers);
    centers = cell2mat(c);
    
elseif Nc1 ~= 2
    error('centers must be a 1xN or 2xN vector!');
end
clear Nc1

cc = ceil(box_size/2);

centers(1,:) = centers(1,:) - cc(1);
centers(2,:) = centers(2,:) - cc(2);
N_centers = size(centers,2);

end_centers = zeros(2,N_centers);
end_centers(1,:) = centers(1,:) + box_size(1)-1;
end_centers(2,:) = centers(2,:) + box_size(2)-1;

for i_c = 1:size(centers,2)
    Ibw(max(centers(1,i_c),1):min(end_centers(1,i_c),size_I(1)),max(centers(2,i_c),1):min(end_centers(2,i_c),size_I(2))) = false;
end


% This seemed more elegant for dealing with boxes out of the image but it
% turned out to be much slower. Computing the canonical window and
% extracting indices like that is still faster if the box is known to be
% fully within the image
% if ~exist('canonical_row_col','var')
%     canonical_center = ceil(box_size/2);
%     canonical_row_col = cell(1,2);
%     canonical_row_col{1} = reshape(repmat((1:box_size(1)) - canonical_center(1),box_size(2),1),[],1)';
%     canonical_row_col{2} = repmat((1:box_size(2)) - canonical_center(2),1,box_size(1));
% end

% % Computing the result:
% res_row = bsxfun(@plus,centers(1,:)',canonical_row_col{1});
% res_col = bsxfun(@plus,centers(2,:)',canonical_row_col{2});
% 
% % Checking sizes:
% valid_ind = res_row > 0 & res_row <= size_I(1) & res_col > 0 & res_col <= size_I(2);
% idx = sub2ind(size_I,res_row(valid_ind),res_col(valid_ind));
% Ibw(idx) = false;
