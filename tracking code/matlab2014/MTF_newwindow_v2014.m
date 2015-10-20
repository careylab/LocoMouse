function [tracks_final,tracks_tail,ONG,tracks_bottom,M,tracks_top,M_top,nong_vect,Unary,Pairwise,xvel,occluded_distance] = MTF_newwindow(data, model)
% MTF   Tracks a set of predefined (mouse) features over a given video.
%
% INPUT:
% data: The data structure
%
% model: Models for detection
%
% This method relies on predefined (or trained) detectors that are specific
% for each feature to be tracked. Features can be either point features
% (e.g. paws, snout) or area features (e.g. tail).
%
% Images are filtered using the provided detectors and post-processing is
% done to extract the candidates:
% For point features, non-maxima suppression is performed a set of discrete
% image positions are kept as candidates.
%
% For area features, contiguous regions that have a detection score above a
% certain threshold are kept as candidates.
%
% Point features are then tracked over the whole sequence in a batch
% (offline) approach using:
%
% Efficient Second Order Multi-Target Tracking with Exclusion Constraints
% Russel C., Setti F., Agapito L.
% British Machine Vision Conference (BMVC 2011), Dundee, UK, 2011.
%
% Code for that tracker is available here:
% http://www.eecs.qmul.ac.uk/~chrisr/tracking.tar.gz.
%
% Area features are not tracked over time. the highest scoring area from
% each image is kept.
%
% Author: Joao Fayad (joaofayad@gmail.com)
% Created: 15/10/2013
% Modified: 06/01/2014
% Modified: 14/10/2015: Cleaning up for release.
%
% The convention about views should be:
% View 1: Bottom.
% View 2: Side.
%
% If it sounds counterintuitive remember that the bottom view is more
% reliable than the side view.
%
% 
% Copyright (C) 2015  Joao Fayad
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.


% Assuming a list of images without background and undistorted:

% Parsing the data structure to something that is better used in parfor:
data_list = data.sequence;
mirror_line = data.mirror_line;
N_images = size(data_list,1);

if N_images == 0
    warning('No images in current dataset!');
    tracks_final = NaN(3,N_images);
    tracks_tail = NaN(3,N_images);
    return;
end

% Area detectors (at this point, just the tail):
w_tail = model.tail.w;
rho_tail = model.tail.rho;
box_size_tail = model.tail.box_size;
detection_threshold_tail = model.tail.detection_threshold;
N_tail_points = 15; % Should be defined somewhere else...

% Point detectors:
w_point = {model.paw.w,model.snout.w};
rho_point = cat(2,model.paw.rho,model.snout.rho);
detection_threshold_point = [model.paw.detection_threshold, model.snout.detection_threshold];
box_size_point = cat(1,model.paw.box_size,model.snout.box_size);


Npt = length(w_point);

% Structures that remain after the parfor loop:
tracks_tail = zeros(3,N_tail_points,N_images);
tracks_joint = cell(Npt,N_images);
tracks_bottom = cell(Npt,N_images);
C = cell(Npt,N_images);

Unary = cell(Npt,N_images);
Pairwise = cell(Npt,N_images-1);
% Tertiary = cell(Npt,N_images-2);

%% Estimating Mouse size and detecting the tail:
% Since the tail is detected independently from the rest of the mouse (for
% now at least) we parse every image and estimate the tail and animal size.
% This allows detecting the mouse with a fixed size box which increases
% robustness.
BB = NaN(2,3,N_images);
Centroids = NaN(2,2,N_images);
for i_images = 1:N_images
    % for i_images = 1:N_images
    I = imreadGrayscaleBrightnessAdjustment(data_list(i_images,:));
    % Estimating bounding box size:
    Ibw = im2bw(medfilt2(I,[10 10]),0.01); % filtering to remove noise.
    I_cell = cell(1,2);
    [I_cell{[2 1]}] = splitImage(Ibw,mirror_line);
    
    bbi = NaN(4,2);
    
    cc = NaN(2,2);
    
    for i_v = 1:2
        
        % Selecting the largest object on the image (and hoping the mouse
        % does not have holes in it...):
        CC = bwconncomp(I_cell{i_v});
        
        [~,largest_object] = max(cellfun(@(x)(length(x)),CC.PixelIdxList));
        % Creating the new Ibw:
        Ibw = false(size(I_cell{i_v}));Ibw(CC.PixelIdxList{largest_object}) = true;
        
        x_sum = sum(Ibw,1);
        y_sum = sum(Ibw,2);
        anybw1 = x_sum > 0;
        anybw2 = y_sum > 0;
        
        cc(:,i_v) = [((1:size(Ibw,2))*x_sum')/sum(x_sum);((1:size(Ibw,1))* y_sum)/sum(y_sum)];
        
        xi = find(anybw1,1,'first');
        xe = find(anybw1,1,'last');
        
        yi = find(anybw2,1,'first');
        ye = find(anybw2,1,'last');
        
        % Saving as BR corner and box size. Note that while the box size
        % is a global measure, the BR corner is defined per image. We chose the
        % right side as mice walk from left to right and so the right side of
        % the box is always the real location (it is not occluded). Remember
        % that images are reversed if animals walk right to left. We choose BR
        % instead of TR so that mapping from TR corner to the others is always
        % done by subtracting the hight and/or width of the box.
        if all(~cellfun(@isempty,{xi,xe,yi,ye}))
            bbi(:,i_v) = [xe;xe-xi;ye;ye-yi];
        end
    end
    Centroids(:,:,i_images) = cc;
    BB(:,:,i_images) = [max(bbi(1:2,:),[],2) bbi(3:4,:)];
end

C2 = reshape(Centroids,4,[]);C2(1,:) = mean(C2([1 3],:),1);C2(3,:) = [];

BOX = median(squeeze(BB(2,:,:)),2);
STD = std(squeeze(BB(2,:,:)),0,2);
BBB = BOX + 3*STD;clear BOX STD % Could take the max, but this might avoid crazy situations...
box_dim = round(min(max(squeeze(BB(2,:,:)),[],2),BBB));
BB = medfilt1(squeeze(BB(1,:,:))',3)';
xvel = diff(squeeze(BB(1,:)));

% Parameters for the tracker:
% Occluded points:
occluded_distance = 15; % Maximum allowed displacement (in pixels).
grid_spacing = 20;
tail_x_threshold = round(box_dim(1)/2);

godo2 = round(grid_spacing/2);

[X,Y] = meshgrid(godo2:grid_spacing:box_dim(2)-godo2,godo2:grid_spacing:round(0.75*box_dim(1)));% Grid of occlusion points.
ONG = repmat([Y(:)';X(:)'],1,1,N_images);

clear X Y
Nong = size(ONG,2);
% alpha_acc = 0;
alpha_vel = 1E-1; % Weight of the velocity term.
Ncandidates = zeros(Npt,N_images);

% Looping over all the images
parfor i_images = 1:N_images
    % Auxiliary variables for the parfor loop:
    I_cell = cell(1,2);
    I = imreadGrayscaleBrightnessAdjustment(data_list(i_images,:));
    [I_cell{[2 1]}] = splitImage(I,mirror_line);
    x_cut = max(BB(1,i_images) - box_dim(1) + 1,1):BB(1,i_images);
    l_x_cut = length(x_cut);
    
    if ~isempty(I_cell{1})
        I_cell{1} = I_cell{1}(max(1,BB(2,i_images) - box_dim(2) + 1):BB(2,i_images),x_cut);
    end
    if ~isempty(I_cell{2})
        I_cell{2} = I_cell{2}(max(1,BB(3,i_images) - box_dim(3) + 1):BB(3,i_images),x_cut);
    end
    
    temp = ONG(:,:,i_images);
    temp(1,:) = BB(1,i_images) - temp(1,:);
    temp(2,:) = mirror_line + BB(2,i_images) - temp(2,:);
    ONG(:,:,i_images) = temp;
    
    OFFSET = BB(:,i_images) - box_dim;
    
    if l_x_cut > tail_x_threshold && ~all(cellfun(@isempty,w_tail))
        % First approximation in searching for the tail only if 50% of the
        % animal is visible:
        x_tail_cut = 1:length(x_cut)-round(0.8*tail_x_threshold);
        
        % Croping the mouse:
        [tracks_tail(:,:,i_images),tm] = tailDetectionFilter3(cellfun(@(x)(x(:,x_tail_cut)),I_cell,'un',0),mirror_line,box_size_tail,w_tail,rho_tail,detection_threshold_tail);
        Tmask = cellfun(@(x,y)([x true(size(x,1),size(y,2)-size(x,2))]),tm,I_cell,'un',0);
        tracks_tail(:,:,i_images) = bsxfun(@plus,tracks_tail(:,:,i_images),bsxfun(@max,OFFSET,zeros(3,1)));
    else
        Tmask{1} = true(size(I_cell{1}));
        Tmask{2} = true(size(I_cell{2}));
        tracks_tail(:,:,i_images) = NaN(3,N_tail_points);
    end
    
    %% Run point detectors:
    % Initializing auxiliary structures:
    scores = cell(1,2);
    i = cell(1,2);
    j = cell(1,2);
    D2 = cell(1,2);
    
    for i_point = 1:Npt
        %% Run filter for each view:
        for i_views = 1:2
            if ~isempty(I_cell{i_views})
                Cmat = conv2(double(I_cell{i_views}), w_point{i_point}{i_views}, 'same') - rho_point{i_views,i_point};
                scores{i_views} = Cmat .* double(im2bw(sc(I_cell{i_views}),0.1)); Cmat(Cmat<0) = 0;
                scores{i_views}(~Tmask{i_views}) = detection_threshold_point(i_views, i_point)-1;
                detections = scores{i_views} > detection_threshold_point(i_views, i_point);
            else
                detections = false;
            end
            if any(detections(:))
                scores{i_views} = scores{i_views}(detections);
                ind_temp = 1:numel(I_cell{i_views});
                % Clustering positions: peakClustering results in the same
                % algorithm as maxg on bbNMS.
                [i{i_views}, j{i_views}] = ind2sub(size(I_cell{i_views}),ind_temp(detections)');
                if i_views  == 1
                    [D2{i_views}, scores{i_views}] = nmsMax([i{i_views} j{i_views}],box_size_point{i_point, i_views}, scores{i_views}','center');
                    D2{i_views} = weightedMean(D2{i_views}, Cmat, box_size_point{i_point,i_views});
                    D2{1}(:,2) = D2{1}(:,2) + max(OFFSET(1),0);
                    D2{1}(:,1) = D2{1}(:,1) + mirror_line + max(OFFSET(2),0);
                    tracks_bottom{i_point,i_images} = [D2{1}(:,[2 1])';scores{1}];
                    Ncandidates(i_point,i_images) = size(tracks_bottom{i_point,i_images},2);
                else
                    [D2{i_views}, scores{i_views}] = peakClustering([i{i_views} j{i_views}],box_size_point{i_point, i_views}, scores{i_views}','center','max');
                    D2{2}(:,1) = D2{2}(:,1) + max(OFFSET(3),0);
                    D2{2}(:,2) = D2{2}(:,2) + max(OFFSET(1),0);
                end
                
            else
                D2{i_views} = zeros(0,2);
                scores{i_views} = zeros(1,0);
                if i_views == 1
                    tracks_bottom{i_point,i_images} = zeros(3,0);
                end
            end
        end
        % Look for detection matches between the views:
        if any(cellfun(@isempty,D2))
            tracks_joint{i_point,i_images} = zeros(5,0);
        else
            % Pairing boxes that overlap horizontally for at least T of box width:
            T = 0.7;
            if i_images > 1
                I_vel = double(imread(data_list(i_images,:))) - double(imread(data_list(i_images-1,:)));
                moving = I_vel > 25;
                D21_mov = reshape(getWindowFromImage(moving,D2{1}',round(box_size_point{i_point,1})/2),[],size(D2{1},1));
                D22_mov = reshape(getWindowFromImage(moving,D2{2}',round(box_size_point{i_point,2})/2),[],size(D2{2},1));
                
                D21_mov = sum(D21_mov,1) >= 0.02*prod(box_size_point{i_point,1});
                D22_mov = sum(D22_mov,1) >= 0.05*prod(box_size_point{i_point,2});
            else
                D21_mov = true(1,size(D2{1},1));
                D22_mov = true(1,size(D2{2},1));
            end

            [Pairs, scores_pairs, C{i_point,i_images}] = boxPairingsWithVelConstraint(D2{1}, D2{2}, scores{1}, scores{2}, D21_mov, D22_mov, box_size_point{i_point, i_views}(2), 'bottom', T);
            tracks_joint{i_point, i_images} = [Pairs scores_pairs]';
            tracks_joint{i_point, i_images} = tracks_joint{i_point, i_images}([2 1 3 4 5],:);
        end
        
        %% Computing unary potentials for tracking:
        % Defining weights for tracks:
        
        switch i_point
            case 1
                weights = zeros(4,size(tracks_bottom{i_point,i_images},2));
                if ~isempty(weights)

                    X = OFFSET(1) + [0 box_dim(1)];
                    Y = OFFSET(2)+ mirror_line + [0 box_dim(2)];
                    weights(4,:) = pawQuadrantWeights_Distance(tracks_bottom{i_point,i_images}(1,:),tracks_bottom{i_point,i_images}(2,:),X,Y,0,0,0.6); % Hind Left
                    weights(2,:) = pawQuadrantWeights_Distance(tracks_bottom{i_point,i_images}(1,:),tracks_bottom{i_point,i_images}(2,:),X,Y,0,1,0.6); % Hind Right
                    weights(3,:) = pawQuadrantWeights_Distance(tracks_bottom{i_point,i_images}(1,:),tracks_bottom{i_point,i_images}(2,:),X,Y,1,0,0.6); % Front Left
                    weights(1,:) = pawQuadrantWeights_Distance(tracks_bottom{i_point,i_images}(1,:),tracks_bottom{i_point,i_images}(2,:),X,Y,1,1,0.6); % Front Right
                    weights = weights';
                    Unary{i_point,i_images} = [repmat(tracks_bottom{i_point,i_images}(3,:)',1,4).*weights;zeros(Nong,4)];
                else
                    Unary{i_point,i_images} = zeros(Nong,4);
                end
                
            case 2
                weights = zeros(1,size(tracks_bottom{i_point,i_images},2));
                if ~isempty(weights)
                    X = OFFSET(1) + [0 box_dim(1)];
                    Y = OFFSET(2)+ mirror_line + [0 box_dim(2)];
                    weights(1,:) = pawQuadrantWeights_Distance(tracks_bottom{i_point,i_images}(1,:),tracks_bottom{i_point,i_images}(2,:),X,Y,1,0.5,0.20); % Front Right
                    weights = weights';
                    valid_snout = weights > 0;
                    tracks_bottom{i_point,i_images} = tracks_bottom{i_point,i_images}(:,valid_snout);
                    Unary{i_point,i_images} = [tracks_bottom{i_point,i_images}(3,:)'.*weights(valid_snout,:);zeros(Nong,1)];
                    Ncandidates(i_point,i_images) = sum(valid_snout);
                    
                else
                    Unary{i_point,i_images} = zeros(Nong,1);
                end
        end
    end
end
%% Computing pairwise potentials:
% Due to the way parfor works we cannot access entry i and i+1 of the same
% cell array without huge overhead. Therefore we sacrifice a bit of memory
% and create an auxiliary structure where each index accesses two frames at
% the time.
clear tracks_tail_bottom tracks_tail_top
% parfor overhead only compensates if sequence is relatively large (~400
% frames).
tracks_bottom_aux = tracks_bottom(:,2:end);
data_list_aux = data_list(2:end,:);

parfor i_images = 1:(N_images-1)
    %     I1 = imread(data_list(i_images,:));
    %     I2 = imread(data_list_aux(i_images,:));
    for i_point = 1:2
        Pairwise{i_point, i_images} = computePairwiseCost(tracks_bottom{i_point,i_images}(1:2,:),tracks_bottom_aux{i_point,i_images}(1:2,:),ONG(:,:,i_images),abs(xvel(i_images))+occluded_distance,alpha_vel);
    end
end

clear tracks_bottom_aux tracks_bottom_aux2
%% Tracking:
% Perform bottom view tracking:
M = cell(Npt,1);
for i_point = 1:Npt
    switch i_point
        case 1
            O = combinator(4,4,'p');
            N_order = size(O,1);
            Cost = zeros(1,N_order);
            m = cell(1,N_order);
            U = Unary(1,:);
            P = Pairwise(1,:);

            parfor i_o = 1:N_order
                u = cellfun(@(x)(x(:,O(i_o,:))),U,'un',false);
                m{i_o} = match2nd(u,P,[],Nong,0);
                c = computeMCost(m{i_o},u,P(1,:));
                Cost(i_o) = sum(c(:));
            end
            [~,imax] = max(Cost);
            M{i_point}(O(imax,:),:) = m{imax};
            
            clear U P T Cost m u c imax
            
        case 2
            M{i_point} =  match2nd (Unary(2,:), Pairwise(2,:), [],Nong, 0);
    end
    
end
M = cell2mat(M);

% Building structures for top view tracking:
M_top = zeros(size(M));
Unary_top = cell(Npt,N_images);
Pairwise_top = cell(Npt,N_images-1);
occluded_distance_top = 15;
N_tracks = size(M,1);
nong_vect = mirror_line-round(occluded_distance_top/2):-occluded_distance_top:1;
alpha_vel = 100;
Nong_top = length(nong_vect);

tracks_final = NaN(3,N_tracks,N_images);
tracks_top = cell(N_tracks,N_images);
for i_tracks = 1:N_tracks
    if i_tracks < 5
        i_point = 1;
    elseif i_tracks == 5
        i_point = 2;
    else
        i_point = 6;
    end
    Ncandidates_top = zeros(1,N_images);
    tracks_joint2 = cell(1,N_images);
    
    % Computing unary potentials:
    parfor i_images = 1:N_images
        
        if M(i_tracks,i_images) <= size(tracks_bottom{i_point,i_images},2)
            first = tracks_joint{i_point,i_images}(1,:) == tracks_bottom{i_point,i_images}(1,M(i_tracks,i_images));
            second = tracks_joint{i_point,i_images}(2,:) == tracks_bottom{i_point,i_images}(2,M(i_tracks,i_images));
            candidates = first & second;
            if any(candidates)
                tracks_joint2{i_images} = tracks_joint{i_point,i_images}([1:3 5],candidates);
                Ncandidates_top(i_images) = sum(candidates);
                Unary_top{i_tracks,i_images} = [tracks_joint2{i_images}(end,:) zeros(1,Nong_top)]';
            else
                tracks_joint2{i_images} = zeros(4,0);
                Unary_top{i_tracks,i_images} = zeros(Nong_top,1);
            end
            
        else
            tracks_joint2{i_images} = zeros(4,0);
            Unary_top{i_tracks,i_images} = zeros(Nong_top,1);
        end
    end
    tracks_top(i_tracks,:) = tracks_joint2;
    
    % Computing pairwise potentials: Abusing a bit of memory to gain speed.
    tracks_joint_2_aux = tracks_joint2(2:N_images);
        
    parfor i_images = 1:(N_images-1)
        Pairwise_top{i_tracks, i_images} = computePairwiseCost(tracks_joint2{i_images}(3,:),tracks_joint_2_aux{i_images}(3,:),nong_vect,occluded_distance_top,alpha_vel);
    end
    clear tracks_joint_2_aux x
    
    % Tracking:
    Mpf_top = match2nd (Unary_top(i_tracks,:), Pairwise_top(i_tracks,:), [],Nong_top, 0);
        
    for i_images = 1:N_images
        
        bottom_cond = M(i_tracks,i_images) <= Ncandidates(i_point,i_images) && M(i_tracks,i_images) > 0;
        top_cond = Mpf_top(i_images) <= Ncandidates_top(i_images) && Mpf_top(i_images) > 0;
        if bottom_cond && top_cond
            tracks_final(:,i_tracks,i_images) = [tracks_bottom{i_point,i_images}(1:2,M(i_tracks,i_images));tracks_joint2{i_images}(3,Mpf_top(i_images))];
        elseif bottom_cond
            tracks_final(:,i_tracks,i_images) = [tracks_bottom{i_point,i_images}(1:2,M(i_tracks,i_images));NaN];
        end
        
        % Performing Z NMS: Perform this step only if more features of the same
        % detector are used (e.g. perform 3 times for paws).
        if i_tracks < 4 && top_cond
            temp = tracks_final(:,i_tracks,i_images)';
            kp = findOverlap(temp(:,[1 3]),tracks_joint{i_point,i_images}([1 3],:)',box_size_point{i_point,2},'center',0.5);
            tracks_joint{i_point,i_images}(:,kp) = repmat(zeros(5,1),1,sum(kp));
        end
    end
    M_top(i_tracks,:) = Mpf_top;
end
