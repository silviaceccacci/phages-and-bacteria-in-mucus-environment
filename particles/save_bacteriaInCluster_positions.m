function bact_in_cluster_positions_over_time = save_bacteriaInCluster_positions(clusters, step_idx, ...
    num_bacteria_total, bact_in_cluster_positions_over_time)
% saveBacteriaInClusterPositions_ID
% ---------------------------------------------------------------
% Stores the positions of bacteria that are part of clusters using fixed IDs.
%
% Inputs:
%   clusters                     : struct array with fields:
%                                    - .size
%                                    - .bacteria: struct array with .position (2x1) and .id
%   step_idx                     : current time step
%   num_bacteria_total           : total number of bacteria in simulation
%   bact_in_cluster_positions_over_time : preallocated matrix
%                                    [num_steps x 2*num_bacteria_total]
%
% Output:
%   bact_in_cluster_positions_over_time : updated matrix
%
% Each row stores positions as:
% [x1, y1, x2, y2, ..., xN, yN] where N = num_bacteria_total
% Bacteria not currently in any cluster are NaN.
% ---------------------------------------------------------------

%     row = NaN(1, 2*max_bacteria_per_cluster*max_num_clusters); % pre-fill with NaNs
% 
%     clusters_filtered = clusters([clusters.size] >= 2);
%     num_clusters_k = numel(clusters_filtered);
% 
%     for c = 1:num_clusters_k
%         if c > max_num_clusters
%             break; % do not exceed preallocation
%         end
% 
%         bacteria_in_c = clusters_filtered(c).bacteria;
%         num_bact_in_c = numel(bacteria_in_c);
% 
%         for b = 1:min(num_bact_in_c, max_bacteria_per_cluster)
%             idx = 2*((c-1)*max_bacteria_per_cluster + (b-1)) + 1;
%             row(idx:idx+1) = bacteria_in_c(b).position(:)'; % assign [x, y]
%         end
%     end

     % Pre-fill with NaNs
    row = NaN(1, 2*num_bacteria_total);

    % Loop over clusters with size >= 2
    clusters_filtered = clusters([clusters.size] >= 2);

    for c = 1:numel(clusters_filtered)
        bacteria_in_c = clusters_filtered(c).bacteria;
        for b = 1:length(bacteria_in_c)
            bact_id = bacteria_in_c(b).id;  % fixed ID
            row(2*bact_id-1 : 2*bact_id) = bacteria_in_c(b).position(:)'; % assign [x, y]
        end
    end

    bact_in_cluster_positions_over_time(step_idx,:) = row;
end
