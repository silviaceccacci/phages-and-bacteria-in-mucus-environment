function cluster_positions_over_time = save_cluster_positions(clusters, step_idx, cluster_positions_over_time)

    max_num_clusters = size(cluster_positions_over_time,2)/2;
    row = NaN(1, 2*max_num_clusters);

    if isempty(clusters)
        cluster_positions_over_time(step_idx, :) = row;
        return;
    end

    sizes = [clusters.size];
    mask = sizes >= 2;
    clusters_filtered = clusters(mask);
    %clusters_filtered = clusters([clusters.size]>=2);
    num_clusters_k = numel(clusters_filtered);

    for c = 1:num_clusters_k
        if c > max_num_clusters
            break;  % avoid exceeding preallocation
        end
        row(2*c-1:2*c) = clusters_filtered(c).position(:)';
    end

    cluster_positions_over_time(step_idx, :) = row;

end
