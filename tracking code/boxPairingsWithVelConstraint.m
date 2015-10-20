function [Pairs, scores, C, sel_ind_bottom, sel_ind_top,D] = boxPairingsWithVelConstraint(D2_bottom, D2_top, scores_bottom, scores_top,mov_bottom, mov_top, box_width, mode, T)
% BOXPAIRINGS Pairs boxes that overlap at least T horizontally.
%
% USAGE: [pairs, scores] = boxPairings(D2_bottom, D2_top, box_widht, T)
%
% INPUT:
% D2_bottom: 2D coordinates of bottom view defined as: [i j].
% D2_bottom: 2D coordinates of top view defined as: [i j].
% scores_bottom: a vector with the score of each bottom detection.
% scores_top: a vector with the score of each top detection.
% box_width: Scalar with the width of both boxes in pixels.
% mode: a string with the merging mode:
%   'weighted': Vertical position is the weighted average (with detect
%   scores).
%   'max': Vertical position is the one with the highest detection score.
% T: Overlapping percentage to consider pairing (between 0 and 1).
%
% Note: The score will always be the sum of both detections, as that is
% what is going to be used during the multi-object tracking. If that
% changes, this must be changed accordingly.
%
%%% FIXME: There is no metric on the detection score, which means that
%%% using them as wheights might create biases. To solve this we need to
%%% normalize scores for something measured at training data.

if ~exist('T','var')
    T = 0.75;
end

overlap_T = round((1-T)*box_width);

N_bottom = size(D2_bottom,1);
N_top = size(D2_top,1);

% W_bottom = repmat(D2_bottom(:,2),1,N_top);
% W_top = repmat(D2_top(:,2)',N_bottom,1);

Distances = pdist2(D2_bottom(:,2),D2_top(:,2));
IsMoving = bsxfun(@eq,mov_bottom(:),mov_top(:)');
% abs(reshape(W_bottom - W_top,N_bottom,N_top));
pairings = (Distances <= overlap_T);% & IsMoving ;

% Check which bottom points are paired more than once, and apply the motion
% constrain on them:
multiple_options = sum(pairings,1)>1;
pairings(:,multiple_options) = pairings(:,multiple_options) & IsMoving(:,multiple_options);

M = reshape(1:numel(pairings),size(pairings));
[sel_ind_bottom,sel_ind_top] = ind2sub(size(pairings),M(pairings));

if any(pairings(:))
    %% Process point pairs:
    lin_indices = reshape(1:numel(pairings),N_bottom,N_top);
    N_pairs = sum(pairings(:));
    
    [i_bottom, i_top] = ind2sub(size(pairings),lin_indices(pairings));
    i_bottom = i_bottom';
    i_top = i_top';
    % correspondences = [i_bottom i_top];
    
    Pairs = zeros(N_pairs,3);
    
    Pairs(:,1) = D2_bottom(i_bottom,1);
    Pairs(:,3) = D2_top(i_top,1);
    % The vertical score is lowered in proportion to the distance that it is
    % moved.
    
    D = ((-Distances(pairings)+overlap_T)/overlap_T);
%     sD = size(D);
%     [~,sD] = max(sD);
%     if sD == 2
%         D = D';
%     end
    sb = scores_bottom(i_bottom);
    st = scores_top(i_top);
    D = D(:);
    scores = [sb(:)  st(:).* D];
    
    % correspondences(correspondences) = 1:size(Pairs,1);
    % [X,Y] = find(pairings);
    % ind = sub2ind(size(pairings),X,Y);
    % pairings(ind) = 1:length(ind);
    P = double(pairings);
    N_pairs = sum(pairings(:));
    P(pairings) = 1:N_pairs;
    
    a = sum(pairings,1) >1;
    C = logical(eye(N_pairs));
    ind = 1:size(pairings,2);
    ind = ind(a);
    for i_c = ind
        who_overlaps = P(pairings(:,i_c),i_c);
        for j_c = 1:length(who_overlaps)
            C(who_overlaps(j_c),who_overlaps(j_c+1:end)) = true;
        end
    end
    
    C = C | C';
    
    
    switch mode
        case 'weighted'
            Pairs(:,2) = round((scores(:,1).* D2_bottom(i_bottom,2) + scores(:,2).*D2_top(i_top,2))./(sum(scores,2)));
            
        case 'maximum'
            Possibilities = [D2_bottom(i_bottom,2) D2_top(i_top,2)];
            
            [~,ind_scores] = max(scores,[],2);
            poss_index = sub2ind(size(Possibilities),1:N_pairs,ind_scores');
            Pairs(:,2) = Possibilities(:,poss_index);
        case 'bottom'
            % For some reason the sizes don't match if we have empty
            % matrices. There is probably another way to do this but this is a
            % quick workaround.
            Pairs(:,2) = reshape(D2_bottom(i_bottom,2),[],1);
            
        otherwise
            error('Unknown option: %s',mode);
    end
else
    Pairs = zeros(0,3);
    scores = zeros(0,2);
    C = false(1,0);
end