function plot_trajectories_2D(xP_over_time, xB_over_time, xC_over_time, num_phages, num_bacteria, num_clusters, micron, Omega_x, Omega_y)
    % Function to plot the 2D trajectories of phages and bacteria over time
    figure;
    hold on;

    % Plot the 2D trajectories of the phages
    for i = 1:num_phages
        %plot(xP_over_time(:,(i-1)*2 + 1) / micron, xP_over_time(:,i*2) / micron, 'bo', 'DisplayName', ['Phage ' num2str(i)], 'LineWidth', 1);
        plot(xP_over_time(:,(i-1)*2 + 1) / micron, xP_over_time(:,i*2) / micron, 'bo');
    end

    % Plot the 2D trajectories of the bacteria
    for j = 1:num_bacteria
        %plot(xB_over_time(:,(j-1)*2 + 1) / micron, xB_over_time(:,j*2) / micron, 'rx', 'DisplayName', ['Bacteria ' num2str(j)], 'LineWidth', 2);
        plot(xB_over_time(:,(j-1)*2 + 1) / micron, xB_over_time(:,j*2) / micron, 'rx');
    end
    
    % Plot the 2D trajectories of the clusters
    for j = 1:num_clusters
        %plot(xB_over_time(:,(j-1)*2 + 1) / micron, xB_over_time(:,j*2) / micron, 'rx', 'DisplayName', ['Bacteria ' num2str(j)], 'LineWidth', 2);
        plot(xC_over_time(:,(j-1)*2 + 1) / micron, xC_over_time(:,j*2) / micron, 'k+');
    end

    % Formatting the plot
    xlabel('$\Omega_x \ (\mu m)$', 'Interpreter','LaTeX','FontSize',16);
    ylabel('$\Omega_y \ (\mu m)$', 'Interpreter','LaTeX','FontSize',16);
    title('2D Trajectories of Phages, Bacteria and Clusters', 'Interpreter','LaTeX','FontSize',16);
    xlim([0 Omega_x/micron]);
    ylim([0 Omega_y/micron]);
    %legend show;
    set(gca,'FontSize',16);
    grid on;
    hold off;

    saveas(gcf, '2D_trajectories' ,'epsc'); 

end
