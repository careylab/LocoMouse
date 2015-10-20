function [Tail,T_mask,tail_length] = tailDetectionFilter3(Ii, mirror_line, box_size, w, rho, detection_thresholds, tail_length)
% TAILDETECTION detects the tail based on thresholding.
%
% INPUT:
% I: I can be a grayscale image or a 1x2 cell. If I is a cell, it must have
% the bottom image on I{1} and the top image on I{2};
% mirror_line: Height at which to split the image in half if I is a
% grayscale image.
% box_size: A 1x2 cell with the box size for the tail on each view (1 -
% bottom, 2 - top).
% w: a 1x2 cell with the filter for the tail on each view (1 - bottom, 2 -
% top).
% rho: a 1x2 cell with the bias value for the tail filte on each view (1 -
% bottom, 2 - top).
%
% OUTPUT:
% Tail:
% T_mask: A mask to remove the tail pixels (and respective boxes) from the
% image.

switch class(Ii)
    case 'uint8'
        I_cell = cell(1,2);
        [I_cell{2:-1:1}] = splitImage(Ii,mirror_line);
    case 'cell'
        I_cell = Ii;
    otherwise
        error('Supported classes for I are uint8 or cell.');
end
clear I

Tail = cell(1,2);
T_mask = cell(1,2);

N = 15;

IB = (conv2(double(I_cell{1}),w{1},'same')-rho{1}); %.* double(im2bw(sc(I_cell{1}),0.1));
IT = (conv2(double(I_cell{2}),w{2},'same')-rho{2}); %.* double(im2bw(sc(I_cell{2}),0.1));

% BB = IB;BB(BB<=detection_thresholds(1)) = 0;BB = logical(BB);
% TT = IT;TT(TT<=detection_thresholds(1)) = 0;TT = logical(TT);
BB = IB > detection_thresholds(1);
TT = IT > detection_thresholds(2);

J = repmat(TT,[1 1 size(IB,1)]) & permute(repmat(BB,[1 1 size(IT,1)]),[3 2 1]);
CC = bwconncomp(J);
if CC.NumObjects > 0
    [~,idmax] = max(cellfun(@(x)(length(x)),CC.PixelIdxList));
    [Z,X,Y] = ind2sub(size(J),CC.PixelIdxList{idmax});
    % Blacking out tail regions:
    ind_bottom = unique(sub2ind(size(I_cell{1}),Y(:),X(:)));
    ind_top = unique(sub2ind(size(I_cell{2}),Z(:),X(:)));
    
    BW = cell(1,2);
    BW{1} = false(size(I_cell{1})); BW{1}(ind_bottom) = true;
    BW{2} = false(size(I_cell{2})); BW{2}(ind_top) = true;
    
    T_mask{1} = true(size(I_cell{1})); T_mask{1} = blackOutImageRegions(T_mask{1},ind_bottom',box_size{1});
    T_mask{2} = true(size(I_cell{2})); T_mask{2}= blackOutImageRegions(T_mask{2},ind_top',box_size{2});
    st = min(X(:));
    en = max(X(:));
    
    % If no length was provided, use the current and divide into N parts:
    if ~exist('tail_length','var')
        tail_length = en-st+1;
    end
    steps = en - (0:N).*round((tail_length)/N);
    steps = steps(steps>0);
%     steps = steps(end:-1:1);
    Ntail = length(steps)-1;
    tail = NaN(2,N);
    
    for i_views = 1:2
        for i_tail = 1:Ntail
            x = sum(BW{i_views}(:,steps(i_tail+1):steps(i_tail)),1);
            y = sum(BW{i_views}(:,steps(i_tail+1):steps(i_tail)),2);
            c = [x*(1:length(x))'/sum(x(:));(1:length(y))*y/sum(y(:))];
            tail(:,i_tail) = [steps(i_tail);0] + round(c);
        end
        
        if i_views == 1
            tail(2,:) = tail(2,:) + mirror_line;
        end
        Tail{i_views} = tail;
    end
else
    Tail = {NaN(2,N),NaN(2,N)};
    T_mask{1} = true(size(I_cell{1}));
    T_mask{2} = true(size(I_cell{2}));
end

T = Tail;
clear Tail
Tail = zeros(3,N);
Tail(1,:) = T{1}(1,:);
Tail(2,:) = T{1}(2,:);
Tail(3,:) = T{2}(2,:);
