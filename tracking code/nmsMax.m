function [D2_nms, scores,kp] = nmsMax(D2_coordinates,box_size, scores, coordinates, overlap)
% NMSMAX_VOL performs non-maxima suppresion in a volume using the Pascal
% criterium of intersection_area/union_area > overlap -> suppression.
%
% USAGE: D2_nms = nmsMax_vol(D2_coordinates, box_size, scores, overlap)
%
% INPUT:
% DD2_coordinates: 2D coordinates defined as: [i j]
% box_size: [h w]
% overlap: decision boundary for supression (between 0 and 1).
% coordinates: string specifying if coordinates are defined at the center
% of the box or at the top left corner (default: center);

% for each i suppress all j st j>i and area-overlap>overlap

if isempty(D2_coordinates)
D2_nms = [];
scores = [];
return;
end

[scores,ord] = sort(scores,'descend'); D2_coordinates = D2_coordinates(ord,:);

if ~exist('coordinates','var')
    coordinates = 'center';
elseif isempty(coordinates)
    coordinates = 'center';
end

if ~any(strcmpi(coordinates,{'center','tlcorner'}))
    error('coordinates must be either ''center'' or ''tlcorner''.');
end

if ~exist('overlap','var')
    overlap = 0.5;
end

n = size(D2_coordinates,1);
kp = true(1,n);
area = prod(box_size);

% Converting center coordinates to box limits:
box_size = repmat(box_size,n,1);
if strcmpi(coordinates,'center')
    DD2_coordinates_s = D2_coordinates - floor(box_size./2);
else
    DD2_coordinates_s = D2_coordinates;
end

DD2_coordinates_e = DD2_coordinates_s + box_size;
clear box_size

for i=1:n
    for j=(i+1):n
        if ~kp(j)
            continue;
        end
        % Overlap:
        iw = min(DD2_coordinates_e(i,1),DD2_coordinates_e(j,1)) - max(DD2_coordinates_s(i,1),DD2_coordinates_s(j,1));
        ih = min(DD2_coordinates_e(i,2),DD2_coordinates_e(j,2)) - max(DD2_coordinates_s(i,2),DD2_coordinates_s(j,2)); 
        if(ih<=0 || iw <=0), continue; end
        
        o = iw * ih;
        u = 2*area - o;
        decision = o/u;
        
        if(decision > overlap)
            kp(j) = false;
        end
    end
end
D2_nms = D2_coordinates(kp,:);
scores = scores(kp);
end