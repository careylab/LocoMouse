function cost = computeMCost(M,Unary,Pairwise)
% COMPUTEMCOST Computes the cost of a transition matrix M given a set of
% Unary and Pairwise terms.

[N_tracks, N_images] = size(M);

cost = zeros(N_tracks,N_images);

for i_images = 1:N_images
    
    for i_track = 1:N_tracks
        
        % Computing Unary term:
        cost(i_track, i_images) = Unary{i_images}(M(i_track,i_images),i_track);
        
        if i_images < N_images
            % Computing Pairwise term:
            cost(i_track, i_images) = cost(i_track, i_images) + Pairwise{i_images}(M(i_track,i_images+1),M(i_track, i_images));
        end
        
    end
    
end
 
