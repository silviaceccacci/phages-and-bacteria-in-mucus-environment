function bact_in_cluster_positions_over_time = save_bacteriaInCluster_positions(clusters, step_idx, ...
    num_bacteria_total, bact_in_cluster_positions_over_time)

    row = NaN(1, 2*num_bacteria_total);

    if isempty(clusters)
        bact_in_cluster_positions_over_time(step_idx, :) = row;
        return;
    end

    sizes = [clusters.size];
    mask = sizes >= 2;
    clusters_filtered = clusters(mask);
    %clusters_filtered = clusters([clusters.size] >= 2);

    for c = 1:numel(clusters_filtered)
        bacteria_in_c = clusters_filtered(c).bacteria;
        for b = 1:length(bacteria_in_c)
            bact_id = bacteria_in_c(b).id;  % fixed ID
            row(2*bact_id-1 : 2*bact_id) = bacteria_in_c(b).position(:)'; % assign [x, y]
        end
    end

    bact_in_cluster_positions_over_time(step_idx,:) = row;

end
