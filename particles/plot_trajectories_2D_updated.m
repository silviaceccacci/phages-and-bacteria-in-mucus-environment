function plot_trajectories_2D_updated(phage_positions_over_time, bact_positions_over_time, cluster_positions_over_time, ...
                              num_phages, num_bacteria, max_num_clusters, micron, Omega_x, Omega_y)
% plot_trajectories_2D
% Plots the 2D trajectories of phages, bacteria, and clusters
%
% Inputs:
%   phage_positions_over_time   : [num_steps x 2*num_phages]
%   bact_positions_over_time    : [num_steps x 2*num_bacteria]
%   cluster_positions_over_time : [num_steps x 2*max_num_clusters] (NaNs allowed)
%   num_phages                  : total number of phages
%   num_bacteria                : total number of bacteria
%   max_num_clusters            : maximum number of clusters
%   micron                      : length scale
%   Omega_x, Omega_y            : domain size

    figure('Color','w'); hold on;

    % ---- Plot bacteria trajectories ----
    for i = 1:num_bacteria
        x = bact_positions_over_time(:, 2*i - 1) / micron;
        y = bact_positions_over_time(:, 2*i) / micron;
        plot(x, y, '-', 'Color', [0.2 0.6 1], 'LineWidth', 1);  % light blue solid line
    end

    % ---- Plot phage trajectories ----
    for i = 1:num_phages
        x = phage_positions_over_time(:, 2*i - 1) / micron;
        y = phage_positions_over_time(:, 2*i) / micron;
        plot(x, y, '-', 'Color', [1 0.4 0.4], 'LineWidth', 1);  % light red solid line
    end

    % ---- Plot cluster COM trajectories ----
    for i = 1:max_num_clusters
        x = cluster_positions_over_time(:, 2*i - 1) / micron;
        y = cluster_positions_over_time(:, 2*i) / micron;
        if all(isnan(x)) && all(isnan(y))
            continue; % skip empty cluster slot
        end
        plot(x, y, 'k--', 'LineWidth', 1.5);  % black dashed line
    end

    % ---- Format plot ----
    xlabel('$x \ (\mu m)$', 'Interpreter','LaTeX','FontSize',16);
    ylabel('$y \ (\mu m)$', 'Interpreter','LaTeX','FontSize',16);
    title('2D Trajectories of Phages, Bacteria and Clusters', 'Interpreter','LaTeX','FontSize',16);
    xlim([0 Omega_x/micron]);
    ylim([0 Omega_y/micron]);
    axis equal;
    grid on;
    set(gca,'FontSize',16);

    legend({'Bacteria','Phages','Clusters (COM)'}, 'Location','best', 'Interpreter','latex');

    hold off;

    % ---- Save figure ----
    %saveas(gcf, 'plot_2D_trajectories','epsc');
end
