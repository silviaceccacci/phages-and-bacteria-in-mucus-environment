function final_cluster_configuration(bacteria, threshold, domain, x_min, y_min, x_max, y_max)
% Forms clusters at final time, prints their info, and plots them

    % 1. Form final clusters
    [final_clusters, bacteria] = Cluster.form_clusters(bacteria, threshold, domain, x_min, y_min, x_max, y_max);

    % 2. Print cluster info
    fprintf('\nFormed %d clusters at final time step\n', length(final_clusters));
    for i = 1:length(final_clusters)
        fprintf('Cluster %d: %d bacteria\n', i, final_clusters(i).size);
        ids = arrayfun(@(b) b.id, final_clusters(i).bacteria);
        fprintf('   Bacteria IDs: %s\n', mat2str(ids));
    end

    % 3. Compute cluster_ids and positionsB_InClu
    num_bacteria = length(bacteria);
    cluster_ids = zeros(num_bacteria, 1);            % bacterium ID → cluster index
    positionsB_InClu = zeros(num_bacteria, 2);       % bacterium ID → position

    for i = 1:length(final_clusters)
        cluster_bacteria = final_clusters(i).bacteria;
        for j = 1:length(cluster_bacteria)
            b_id = cluster_bacteria(j).id;
            cluster_ids(b_id) = i;
            positionsB_InClu(b_id, :) = cluster_bacteria(j).position;
        end
    end

    % Remove unassigned (non-clustered) entries
    valid = cluster_ids > 0;
    cluster_ids = cluster_ids(valid);
    positionsB_InClu = positionsB_InClu(valid, :);

    % 4. Compute center of mass for each cluster
    num_clusters = length(final_clusters);
    comC = zeros(num_clusters, 2);
    for i = 1:num_clusters
        final_clusters(i) = final_clusters(i).update_center_of_mass_and_group_velocity();
        comC(i, :) = final_clusters(i).center_of_mass;
    end

    % 5. Plot final cluster configuration
    colors = lines(num_clusters);  % One color per cluster

    figure;
    hold on;

    for i = 1:num_clusters
        % Get bacteria belonging to cluster i
        idx = find(cluster_ids == i);
        pos = positionsB_InClu(idx, :);

        % Plot bacteria
        scatter(pos(:,1), pos(:,2), 60, 'filled', ...
            'MarkerFaceColor', colors(i,:), ...
            'DisplayName', sprintf('Cluster %d', i));

        % Label bacteria
        for j = 1:length(idx)
            text(pos(j,1) + 0.1, pos(j,2), ...
                sprintf('B%d', idx(j)), ...
                'FontSize', 10, 'Color', colors(i,:));
        end

        % Plot center of mass
        scatter(comC(i,1), comC(i,2), 100, 'p', ...
            'MarkerEdgeColor', 'k', ...
            'MarkerFaceColor', colors(i,:), ...
            'LineWidth', 1.5);
    end

    xlabel('x-position','Interpreter','LaTeX','FontSize',14);
    ylabel('y-position','Interpreter','LaTeX','FontSize',14);
    title('Clusters at Final Time','Interpreter','LaTeX','FontSize',16);
    legend('show');
    axis equal;
    set(gca, 'FontSize', 14);
    hold off;

    saveas(gcf, 'plot_final_clusters', 'epsc');

end
