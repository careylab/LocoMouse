function [D2_pc, scores] = peakClustering(D2_coordinates,box_size, scores, coordinates, mode, overlap)
% PEAKCLUSTERING Clusters detections that overlap. The new center is their
% weighted mean by the detection scores. The new score is the median
% detection score.
%
% USAGE: [D2_pc, scores] = peakClustering(D2_coordinates, box_size, scores, coordinates, mode, overlap)
%
% INPUT:
% D2_coordinates: 2D coordinates defined as: [i j].
% box_size: [h w].
% scores: a vector with the score of each detection.
% overlap: decision boundary for supression (between 0 and 1).
% coordinates: a string specifying if coordinates are given in box center
% or box top left corner notation (default: 'center').
% mode: a string defining the cluster score method as the 'max', 'mean',
% 'median' or 'sum' of the scores of clustered points (default: 'median');

if isempty(D2_coordinates)
    D2_pc = D2_coordinates;
    return;
else
    
    if ~exist('mode','var')
        mode = 'median';
    end
    
    if ~exist('coordinates','var')
        coordinates = 'center';
    end
    
    if ~any(strcmpi(coordinates,{'center','tlcorner'}))
        error('coordinates must be either ''center'' or ''tlcorner''.');
    end
    
    if ~any(strcmpi(mode,{'median','max','mean','sum'}))
        error('mode must be either ''median'', ''max'' or ''mean''.');
    end
    
    if ~exist('overlap','var')
        overlap = 0.5;
    end
    
    % Sort detections according to score:
    [scores,ord] = sort(scores,'descend'); D2_coordinates = D2_coordinates(ord,:);
    
    n = size(D2_coordinates,1);
    kp = true(1,n);
    area = prod(box_size);
    
    % Converting center coordinates to box limits:
    box_size = repmat(box_size,n,1);
    if strcmpi(coordinates,'center')
        D2_coordinates_start = D2_coordinates - floor(box_size./2);
    else
        D2_coordinates_start = D2_coordinates;
    end
    D2_coordinates_end = D2_coordinates_start + box_size;
    
    % load A.mat
    % c = 0;
    for i=1:n
        if ~kp(i)
            continue;
        end
        current_cluster = false(1,n);
        current_cluster(i) = true;
        for j=(i+1):n
            if ~kp(j)
                continue;
            end
            
            % Overlap:
            iw = min(D2_coordinates_end(i,1),D2_coordinates_end(j,1)) - max(D2_coordinates_start(i,1),D2_coordinates_start(j,1));
            ih = min(D2_coordinates_end(i,2),D2_coordinates_end(j,2)) - max(D2_coordinates_start(i,2),D2_coordinates_start(j,2));
            
            % Do not overlap:
            if(ih<=0 || iw <=0), continue; end
            
            % Computing the Pascal criterion:
            o = iw * ih;
            u = 2*area - o;
            decision = o/u;
            
            % Decision:
            if(decision > overlap)
                kp(j) = false;
                current_cluster(j) = true;
            end
        end
        
        D2_coordinates(i,:) = round(sum(D2_coordinates(current_cluster,:).*repmat(scores(current_cluster)',1,2),1)./sum(scores(current_cluster)));
        scores(i) = feval(mode,scores(current_cluster));
        %     ind = sub2ind(size(A),D2_coordinates(kp,1),D2_coordinates(kp,2));
        %     fA = true(size(A));
        %     fA(ind) = false;
        %     A(fA) = 0;
        %     c = c+1;
        %     imwrite(sc(A,[0 1],'jet'),sprintf('im_%d05.png',c));
    end
    D2_pc = D2_coordinates(kp,:);
    scores = scores(kp);
    
    % ind = sub2ind(size(A),D2_pc(:,1),D2_pc(:,2));
    % fA = true(size(A));
    % fA(ind) = false;
    % A(fA) = 0;
    % c = c+1;
    % imwrite(sc(A,[0 1],'jet'),sprintf('im_%d05.png',c));
end