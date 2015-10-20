function [Iwindow,accepted_coordinates] = getWindowFromImage(I,box_center,box_size, partial_view)
%GETWINDOWFROMIMAGE     extracts a region of a given image
%
% GETWINDOWFROMIMAGE(I, box_center, box_size) extracts the subregion of I centered
% in box_center and of size box_size.
% I: a 2D or 3D matrix.
% box_center: i and j coordinates in the image I (can be a stack of centers).
% box_size: height and width of the region to extract.
% partial_view: boolean controlling the validity of boxes that go outside
% the image. True (default) accepts partial detections, with zeros on
% pixels outside the image.
%
% OUTPUT:
% Iwindow: a tensor stacking the boxes extracted.
% accepted_coordinates: boolean vector containing true if a window was
% successfully extracted and false otherwise (success depends on
% partial_view).
% 
% Note: If a box is not successfully extracted, in its place in the tensor
% there will still be a box filled with zeros. Selecting the successfull
% images must be done with accepted_coordinates outside of this function.
%
% Example:
%
% I = imread('peppers.png');
% box_center = [50 100];
% box_size = [200 300];
% I2 = getWindowFromImage(I, box_center, box_size);
% imshow(I2);

classI = class(I);
if isempty(box_center);
    eval(['Iwindow = ' classI '(zeros(box_size(1),box_size(2),0));']);
    accepted_coordinates = true(0);
    return;
end

if ~exist('partial_view','var')
    partial_view = false;
end

if size(box_center,1) ~= 2
    error('box_center must be a 2xN matrix with the i and j coordinates of the windows.');
end

if length(box_size) ~= 2
    error('box_size must be a 2-vector with the i and j coordinates of the window box_center');
end

if any(box_center < 1) | any(box_size < 1)
    error('All entries in box_center and box_size must be greater than 1.');
end

Ncenters = size(box_center,2);

if Ncenters > 0
    
    i_add = box_center(1,:) - ceil(box_size(1)/2);
    i_vec = repmat(1:box_size(1),Ncenters,1) + repmat(i_add',1,box_size(1));
    
    j_add = box_center(2,:) - ceil(box_size(2)/2);
    j_vec = repmat(1:box_size(2),Ncenters,1) + repmat(j_add',1,box_size(2));
    
    size_I = size(I);
    
    eval(['Iwindow = ' classI '(zeros([box_size(1),box_size(2),Ncenters]));']);
    
    ind_row = 1:box_size(1);
    ind_col = 1:box_size(2);
    accepted_coordinates = true(1,Ncenters);
    
    for n = 1:Ncenters
        row_map = (i_vec(n,:) > 0 & i_vec(n,:) <= size_I(1));
        col_map = (j_vec(n,:) > 0 & j_vec(n,:) <= size_I(2));
        
        if ~partial_view && ~all(row_map) && ~all(col_map)
            accepted_coordinates(n) = false;
        else
            Iwindow(ind_row(row_map),ind_col(col_map),n) = I(i_vec(n,row_map),j_vec(n,col_map));
        end
    end
else
    eval(['Iwindow = ' classI '(zeros(box_size(1),box_size(2),0));']);
end
