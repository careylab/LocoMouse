function kp = findOverlap(D2_candidate,D2_coordinates,box_size,coordinates,overlap)
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

if isempty(D2_candidate)
    kp = false(1,size(D2_coordinates,2));
    return;
end

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
kp = false(1,n);
area = prod(box_size);

% Converting center coordinates to box limits:
box_size = repmat(box_size,n,1);
if strcmpi(coordinates,'center')
    DD2_coordinates_s = D2_coordinates - floor(box_size./2);
    D2_candidate_s = D2_candidate - floor(box_size(1)/2);
else
    DD2_coordinates_s = D2_coordinates;
    D2_candidate_s = D2_candidate;
end

DD2_coordinates_e = DD2_coordinates_s + box_size;
D2_candidate_e = D2_candidate_s + box_size(1,:);
clear box_size

for j=1:n
    % Overlap:
    iw = min(D2_candidate_e(1),DD2_coordinates_e(j,1)) - max(D2_candidate_s(1),DD2_coordinates_s(j,1));
    ih = min(D2_candidate_e(2),DD2_coordinates_e(j,2)) - max(D2_candidate_s(2),DD2_coordinates_s(j,2));
    if(ih<=0 || iw <=0), continue; end
    
    o = iw * ih;
    u = 2*area - o;
    decision = o/u;
    
    if(decision > overlap)
        kp(j) = true;
    end
end