function tracks_tail= correctTracksFlip2(tracks_tail, image_size, mirror_line)
% CORRECTTRACKSFLIP Coordinate transformation for trials where mice run
% from right to left.
%
% The mouse tracking code is based on the assumption that mice run from
% left to right (due to legacy and implementation reasons). When mice run
% from right to left, the image is flipped before tracking. These
% coordinates must then be corrected so they can be analysed just like the
% trials where mice run from left to right.
%
% The 180deg rotation is equivalent to flipping the X and Y coordinates.
% Since X is already flipped for tracking, We only need to flip Y. Z is not
% affected by the change in camera position.
%
% Additionally, we need to reverse the L/R limb labeling as the algorithm
% labels the closest side as R but in flipped images the correct side is L.
 
[~, N_tracks, ~] = size(tracks_tail);

%% 1) Correcting the Y coordinate
flip_coord = [zeros(1,mirror_line) image_size(1):-1:mirror_line+1];
T = ~isnan(tracks_tail(2,:));
%final_tracks(2,T) = flip_coord(final_tracks(2,T));
tracks_tail(2,T) = flip_coord(tracks_tail(2,T));

%% 2) Switch Left/Right labelling:
%final_tracks = final_tracks(:,[3 4 1 2 N_tracks-4:N_tracks],:);
