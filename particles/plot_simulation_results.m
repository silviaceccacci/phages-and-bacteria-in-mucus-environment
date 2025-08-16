function plot_simulation_results(time_vec, cluster_sizes_over_time, phage_attachement_history, DeltaK)
%PLOT SIMULATION RESULTS Plots clustering and phage attachment metrics over time and saves summary data.
%
% Inputs:
%   time_vec - Time vector [1 x T]
%   cluster_sizes_over_time - Cell array of cluster sizes at each time step {1 x T}
%   phage_attachement_history - Matrix of phage attachments [T x N_bacteria]
%   DeltaK - scalar representing difference in porosity

%% 1. Plot: Mean number of clusters over time
    figure;
    num_clusters_over_time = cellfun(@numel, cluster_sizes_over_time);
    plot(time_vec, num_clusters_over_time, 'LineWidth', 1.5);
    xlabel('Time (s)','Interpreter','LaTeX','FontSize',16);
    ylabel('Number of clusters','Interpreter','LaTeX','FontSize',16);
    title('Cluster count over time','Interpreter','LaTeX','FontSize',16);
    set(gca,'FontSize',16);
    %saveas(gcf, 'plot_clusters_count' ,'epsc'); 
    saveas(gcf, sprintf('plot_clusters_count_DeltaK_%g', DeltaK) ,'epsc'); 

    % Save cluster counts to CSV file [time, cluster_counts]
    cluster_csv_filename = sprintf('DeltaK_%g_cluster_counts.csv', DeltaK);
    csv_data = [time_vec(:), num_clusters_over_time(:)];
    writematrix(csv_data, cluster_csv_filename);

 %% 2. Plot: Phage attachment evolution over time (all bacteria)
    figure;
    plot(time_vec, phage_attachement_history, 'LineWidth',2);
    xlabel('Time (s)','Interpreter','LaTeX','FontSize',16);
    ylabel('Number of attached phages','Interpreter','LaTeX','FontSize',16);
    legend(arrayfun(@(i) sprintf('B%d', i), 1:size(phage_attachement_history,2), 'UniformOutput', false), 'Location', 'eastoutside');
    title('Phage attachments per bacterium over time','Interpreter','LaTeX','FontSize',16);
    set(gca,'FontSize',16);
    saveas(gcf, sprintf('plot_phages_attachments_DeltaK_%g', DeltaK) ,'epsc'); 
    hold off;

%     % 2. Plot: Phage attachment evolution over time
%     non_zero_bacteria = any(phage_attachement_history > 0, 1); % logical array [1 x N_bacteria]
%     filtered_data = phage_attachement_history(:, non_zero_bacteria);
%     filtered_indices = find(non_zero_bacteria); % bacteria indices to use in legend
% 
%     figure;
%     plot(time_vec, filtered_data, 'LineWidth',2);
%     xlabel('Time (s)','Interpreter','LaTeX','FontSize',16);
%     ylabel('Number of attached phages','Interpreter','LaTeX','FontSize',16);
%     legend(arrayfun(@(i) sprintf('B%d', i), filtered_indices, 'UniformOutput', false));
%     title('Phage attachments per bacterium over time','Interpreter','LaTeX','FontSize',16);
%     set(gca,'FontSize',16);
%     saveas(gcf, 'plot_phages_attachments' ,'epsc'); 
%     hold off;

%% 3. Plot: Final phage attachment snapshot for all bacteria (bar plot)
    final_attachment_counts = phage_attachement_history(end, :);  % Last time step

    figure;
    bar(1:length(final_attachment_counts), final_attachment_counts, 'FaceColor', [0.2 0.6 0.8]);
    xlabel('Bacterium ID','Interpreter','LaTeX','FontSize',16);
    ylabel('Number of attached phages','Interpreter','LaTeX','FontSize',16);
    title('Phage attachments at final time step (all bacteria)','Interpreter','LaTeX','FontSize',16);
    xticks(1:length(final_attachment_counts));
    xticklabels(arrayfun(@(i) sprintf('B%d', i), 1:length(final_attachment_counts), 'UniformOutput', false));
    set(gca,'FontSize',16);
    saveas(gcf, sprintf('plot_phages_attachments_final_DeltaK_%g', DeltaK), 'epsc');

%     % 3. Plot: Final phage attachment snapshot (bar plot)
%     final_attachment_counts = phage_attachement_history(end, :);  % Last time step
%     non_zero_final = final_attachment_counts > 0;
%     final_counts_filtered = final_attachment_counts(non_zero_final);
%     final_ids_filtered = find(non_zero_final);
% 
%     figure;
%     bar(final_ids_filtered, final_counts_filtered, 'FaceColor', [0.2 0.6 0.8]);
%     xlabel('Bacterium ID','Interpreter','LaTeX','FontSize',16);
%     ylabel('Number of attached phages','Interpreter','LaTeX','FontSize',16);
%     title('Phage attachments at final time step','Interpreter','LaTeX','FontSize',16);
%     xticks(final_ids_filtered);
%     xticklabels(arrayfun(@(i) sprintf('B%d', i), final_ids_filtered, 'UniformOutput', false));
%     set(gca,'FontSize',16);
%     saveas(gcf, 'plot_phages_attachments_final', 'epsc');
   
    %% 4. Plot: Average number of attached phages per bacterium over time
    avg_phages_per_bacterium = mean(phage_attachement_history, 2); % average over bacteria for each time step

    figure;
    plot(time_vec, avg_phages_per_bacterium, 'LineWidth', 2);
    xlabel('Time (s)','Interpreter','LaTeX','FontSize',16);
    ylabel('$<N_p>$','Interpreter','LaTeX','FontSize',16);
    %ylabel('Average number of attached phages per bacterium','Interpreter','LaTeX','FontSize',16);
    title('Average phage attachment per bacterium over time','Interpreter','LaTeX','FontSize',16);
    set(gca,'FontSize',16);
    saveas(gcf, sprintf('plot_avg_phage_attachment_DeltaK_%g', DeltaK), 'epsc');

    % Save average phage attachment to CSV file [time, avg_phage_attachment]
    avg_phage_csv_filename = sprintf('DeltaK_%g_avg_phage_attachment.csv', DeltaK);
    csv_avg_data = [time_vec(:), avg_phages_per_bacterium(:)];
    writematrix(csv_avg_data, avg_phage_csv_filename);

end
