function pairwise_cost = computePairwiseCost(Xi,Xip1,ONG,occluded_distance,alpha_vel)

dim = size(ONG,1);

% Compute second order costs:
if isempty(Xi)
    Xi = zeros(dim,0);
end

if isempty(Xip1)
    Xip1 = zeros(dim,0);
end

Ni = size(Xi,2);
Nip1 = size(Xip1,2);
Nong = size(ONG,2);

p1 = repmat(Xi,1,Nip1);
p2 = reshape(reshape(repmat(Xip1,Ni,1),[],1),dim,Ni*Nip1);

pair_diff = p1 - p2;
S1 = sqrt(sum(pair_diff.^2,1));
if numel(occluded_distance) == 1
    invalid_transitions = S1 > occluded_distance;
    temp = reshape(S1,Ni,Nip1);
    temp(invalid_transitions) = 0;
    temp(~invalid_transitions) = occluded_distance - temp(~invalid_transitions);
    temp = temp/occluded_distance;
else
    pair_diff = abs(pair_diff);
    invalid_transitions = any(bsxfun(@gt,pair_diff,occluded_distance),1);
    temp = zeros(Ni,Nip1);
    temp2 = bsxfun(@minus,occluded_distance,pair_diff(:,~invalid_transitions));
    temp2 = temp2 ./ repmat(occluded_distance,1,size(temp2,2));
    temp(~invalid_transitions) = sum(temp2,1)/dim;
    clear temp2
end
% clear temp
% Occluded costs:
% The cost of becoming occluded: A point can only transition to its nearest
% occlusion point. The cost of that transition is set here:
idx = knnsearch(ONG',Xi');
temp_occ2 = zeros(Ni,Nong);
idx = sub2ind(size(temp_occ2),1:Ni,idx');
temp_occ2(idx) = 0.01;
temp_occ3 = 0.01 * eye(Nong);

% The cost of becoming visible: A point can only resurface to candidate
% positions assigned to that occlusion point.
idx = knnsearch(ONG',Xip1');
temp_occ = zeros(Nong,Nip1);
idx = sub2ind(size(temp_occ),idx',1:Nip1);
temp_occ(idx) = 0.01;%occluded_distance;

pairwise_cost = sparse(alpha_vel * [temp temp_occ2;temp_occ temp_occ3]');