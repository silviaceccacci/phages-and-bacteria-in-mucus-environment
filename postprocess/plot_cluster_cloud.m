function plot_cluster_cloud(cluster, color, alpha, threshold_distance)
%PLOT_CLUSTER_CLOUD Plots a soft density cloud over a given cluster.
%   cluster : a Cluster object with .bacteria having .position
%   color   : RGB triplet for the patch
%   alpha   : transparency of the patch

    % Extract positions of real bacteria
    pos = cat(1, cluster.bacteria.position);

    % Skip if only one bacterium
    if size(pos, 1) == 1
        return
    end

    % --- Add ghost point if only 2 bacteria to form a plane ---
    if size(pos, 1) == 2
        p1 = pos(1, :); 
        p2 = pos(2, :);
        center = mean([p1; p2]);
        v = p2 - p1;
        v_perp = [-v(2), v(1)];
        v_perp = v_perp / norm(v_perp + eps);  % safe normalize
        ghost_point = center + 1e-3 * v_perp;
        pos_with_ghost = [pos; ghost_point];
    else
        pos_with_ghost = pos;
    end

    % --- Density estimation over a grid ---
    x = pos_with_ghost(:,1); 
    y = pos_with_ghost(:,2);
    margin = threshold_distance/2;
    [X, Y] = meshgrid(linspace(min(x)-margin, max(x)+margin, 100), ...
                      linspace(min(y)-margin, max(y)+margin, 100));

    % --- Perform kernel density estimate ---
    n = size(pos_with_ghost, 1);
    bandwidth = margin/2;  % adjust for tighter/looser cloud
    f = zeros(size(X));
    for i = 1:n
        dx = X - pos_with_ghost(i,1);
        dy = Y - pos_with_ghost(i,2);
        f = f + exp(-(dx.^2 + dy.^2)/(2*bandwidth^2));
    end
    f = f / (2*pi*bandwidth^2 * n);

    % --- Plot a contourf "cloud" where density is significant ---
     contourf(X, Y, f, [max(f(:))/4 max(f(:))/6], ...
             'FaceColor', color, 'FaceAlpha', alpha, 'LineStyle', 'none');
end

