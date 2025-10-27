function cluster_positions_over_time = save_cluster_positions(clusters, step_idx, cluster_positions_over_time)
% saveClusterPositions
% ------------------------------------------------------------------
% Stores the positions (centers of mass) of clusters with size >= 2
% at a given time step into a preallocated matrix.
%
% Usage:
%   cluster_positions_over_time = saveClusterPositions(clusters, step_idx, cluster_positions_over_time)
%
% Inputs:
%   clusters                    : struct array with fields .size and .position (2×1 vector)
%   step_idx                    : current time step index
%   cluster_positions_over_time : preallocated matrix [num_steps, 2 * max_num_clusters]
%
% Output:
%   cluster_positions_over_time : updated matrix with current cluster COMs stored
%
% Each row of cluster_positions_over_time stores:
%   [x1, y1, x2, y2, ..., xM, yM]
% where M is the number of clusters with size ≥ 2 at that step.
% Unused entries remain NaN.
%
% ------------------------------------------------------------------

    % Filter clusters with size >= 2
%     clusters_filtered = clusters([clusters.size] >= 2);
%     num_clusters_k = numel(clusters_filtered);
% 
%     % Pre-fill row with NaNs
%     row = NaN(1, size(cluster_positions_over_time, 2));
% 
%     for c = 1:num_clusters_k
%         pos = clusters_filtered(c).position(:); % 2x1
%         idx = 2*c - 1;
%         row(idx:idx+1) = pos;
%     end

    max_num_clusters = size(cluster_positions_over_time,2)/2;
    row = NaN(1, 2*max_num_clusters);

    clusters_filtered = clusters([clusters.size]>=2);
    num_clusters_k = numel(clusters_filtered);

    for c = 1:num_clusters_k
        if c > max_num_clusters
            break;  % avoid exceeding preallocation
        end
        row(2*c-1:2*c) = clusters_filtered(c).position(:)';
    end

    cluster_positions_over_time(step_idx, :) = row;

%     if num_clusters_k > 0
%         % Extract COM positions into a 2×M matrix
%         xC_yC_matrix = NaN(2, num_clusters_k);
%         for c = 1:num_clusters_k
%             xC_yC_matrix(:, c) = clusters_filtered(c).position(:);
%         end
% 
%         % Flatten as [x1, y1, x2, y2, ...]
%         xC_yC_row_fixed_time = reshape(xC_yC_matrix, [2*num_clusters_k, 1])';
% 
%         % Assign to current time step (preserve previous NaNs beyond current cluster count)
%         cluster_positions_over_time(step_idx, 1:length(xC_yC_row_fixed_time)) = xC_yC_row_fixed_time;
%     end
end
