function [coeff, score, explained, global_mean] = compute_globalClusterPCA(coordC_over_time, cluster_sizes_over_time)
% Perform PCA on the trajectories of cluster centers.
%
% Input:
%   coordC_over_time : [N_steps x 2*max_clusters] matrix where each row
%                      contains [x1 y1 x2 y2 ... xM yM] cluster positions at one time step.
%
% Output:
%   coeff    : Principal component coefficients (directions)
%   score    : Transformed data (projection onto principal components)
%   latent   : Eigenvalues (variance along principal components)
%   explained: Percentage of variance explained by each component


% compute_globalClusterPCA - Performs PCA on the trajectories of all clusters with size >= 2
%
% Inputs:
%   cluster_positions_over_time - cell array of cluster positions over time (length T, each cell is a Cx2 matrix)
%   cluster_sizes_over_time     - cell array of cluster sizes over time (length T, each cell is a Cx1 vector)
%
% Outputs:
%   coeff     - Principal directions (2x2)
%   score     - Transformed coordinates (Nx2)
%   explained - Percentage of variance explained by each PC
%   global_mean - Mean of all points used for PCA

    N_steps = size(coordC_over_time, 1);
    cluster_positions = [];

    for k = 1:N_steps
        sizes_k = cluster_sizes_over_time{k};
        num_clusters_k = length(sizes_k);  % clusters with size >= 2 at step k

        if num_clusters_k == 0
            continue;
        end

        row_k = coordC_over_time(k, 1:(2*num_clusters_k));  % get active cluster coordinates at step k
        positions_k = reshape(row_k, [2, num_clusters_k])';  % (num_clusters_k x 2)
        cluster_positions = [cluster_positions; positions_k];  % accumulate all positions
    end

    if isempty(cluster_positions)
        error('No clusters found with size >= 2 for PCA.');
    end

    % Subtract mean for PCA
    global_mean = mean(cluster_positions, 1);
    centered_positions = cluster_positions - global_mean;

    % Perform PCA
    [coeff, score, latent] = pca(centered_positions);

    % Variance explained
    explained = 100 * latent / sum(latent);

    % Optional: plot result
    figure;
    scatter(cluster_positions(:,1), cluster_positions(:,2), 30, 'b', 'filled');
    hold on;
    quiver(global_mean(1), global_mean(2), coeff(1,1), coeff(2,1), ...
           'r', 'LineWidth', 2, 'MaxHeadSize', 2);
    title(sprintf('Global PCA - PC1 explains %.2f%%', explained(1)));
    xlabel('x [µm]');
    ylabel('y [µm]');
    axis equal;
    grid on;

    saveas(gcf, 'global_pca' ,'epsc'); 

end
