function plot_all_trajectories(phage_positions_over_time, bact_positions_over_time, cluster_positions_over_time, time_vec, outputFolder)
% plotAllTrajectories
% -------------------------------------------------------------------------
% Plots trajectories of phages, bacteria, and clusters COMs from recorded
% position matrices.
%
% Usage:
%   plotAllTrajectories(phage_positions_over_time, bact_positions_over_time, ...
%                       cluster_positions_over_time, time_vec, outputFolder)
%
% Inputs:
%   phage_positions_over_time   : [num_steps × 2*num_phages]
%   bact_positions_over_time    : [num_steps × 2*num_bacteria]
%   cluster_positions_over_time : [num_steps × 2*max_num_clusters] (NaNs allowed)
%   time_vec                    : [num_steps × 1] vector of times
%   outputFolder                : folder to save figure (string)
%
% Output:
%   Saves a figure 'trajectories_all.png' in outputFolder
%
% -------------------------------------------------------------------------

    % Create figure
    figure('Color', 'w'); hold on;

    % ---- Plot bacteria trajectories ----
    num_bacteria = size(bact_positions_over_time, 2) / 2;
    for i = 1:num_bacteria
        x = bact_positions_over_time(:, 2*i - 1);
        y = bact_positions_over_time(:, 2*i);
        plot(x, y, '-', 'Color', [0.2 0.6 1], 'LineWidth', 0.8); % light blue
    end

    % ---- Plot phage trajectories ----
    num_phages = size(phage_positions_over_time, 2) / 2;
    for i = 1:num_phages
        x = phage_positions_over_time(:, 2*i - 1);
        y = phage_positions_over_time(:, 2*i);
        plot(x, y, '-', 'Color', [1 0.4 0.4], 'LineWidth', 0.8); % light red
    end

    % ---- Plot cluster COM trajectories ----
    num_clusters_max = size(cluster_positions_over_time, 2) / 2;
    for i = 1:num_clusters_max
        x = cluster_positions_over_time(:, 2*i - 1);
        y = cluster_positions_over_time(:, 2*i);
        if all(isnan(x)) && all(isnan(y))
            continue; % skip empty cluster slot
        end
        plot(x, y, 'k--', 'LineWidth', 1.5); % black dashed line
    end

    % ---- Format ----
    xlabel('$x$ ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
    ylabel('$y$ ($\mu m$)', 'Interpreter', 'latex', 'FontSize', 14);
    title('Trajectories of Phages, Bacteria, and Cluster COMs', 'Interpreter', 'latex', 'FontSize', 14);
    legend({'Bacteria', 'Phages', 'Clusters (COM)'}, 'Location', 'best', 'Interpreter', 'latex');
    axis equal;
    grid on;
    box on;

    % ---- Save figure ----
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end
    saveas(gcf, fullfile(outputFolder, 'trajectories_all.png'));
    hold off;

    disp(['✅ Saved trajectory plot to ' fullfile(outputFolder, 'trajectories_all.png')]);
end
