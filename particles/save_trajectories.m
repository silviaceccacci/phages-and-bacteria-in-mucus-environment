function save_trajectories(xP_over_time, xB_over_time, xC_over_time, dt, N_steps, num_phages, num_bacteria, num_clusters, micron, outputFolder)

    %% Save Phage Positions to a .dat File
    phage_filename = fullfile(outputFolder, 'trajectoriesPhages.dat'); 
    phage_file = fopen(phage_filename, 'w'); 

    % Write header for phage file
    fprintf(phage_file, 'Time (s) \t');
    for i = 1:num_phages
        fprintf(phage_file, 'P%d X (micron) \t P%d Y (micron) \t', i, i);
        %fprintf(phage_file, 'Phage %d X (micron) \t Phage %d Y (micron) \t', i, i);
    end
    fprintf(phage_file, '\n');

    % Write the phage positions at each time step
    for k = 1:N_steps
        fprintf(phage_file, '%.6f \t', (k-1)*dt); % Time
        for i = 1:num_phages
            % Extract the X and Y positions for phage i at time step k
            x_val = xP_over_time(k, (i-1)*2 + 1) / micron;
            y_val = xP_over_time(k, i*2) / micron;

            %fprintf('Time: %.6f, Phage %d X: %.6f, Y: %.6f\n', (k-1)*dt, i, x_val, y_val);
            
            fprintf(phage_file, '%.6f \t %.6f \t', x_val, y_val);
        end
        fprintf(phage_file, '\n');
    end
    fclose(phage_file); 
    
    %% Save Bacteria Positions to a .dat File
    bacteria_filename = fullfile(outputFolder, 'trajectoriesBacteria.dat'); 
    bacteria_file = fopen(bacteria_filename, 'w'); 

    % Write header for bacteria file
    fprintf(bacteria_file, 'Time (s)\t');
    for j = 1:num_bacteria
        fprintf(bacteria_file, 'B%d X (micron) \t B%d Y (micron) \t', j, j);
        %fprintf(bacteria_file, 'Bacteria %d X (micron) \t Bacteria %d Y (micron) \t', j, j);
    end
    fprintf(bacteria_file, '\n');

    % Write the bacteria positions at each time step
    for k = 1:N_steps
        fprintf(bacteria_file, '%.6f \t', (k-1)*dt); % Time
        for j = 1:num_bacteria
            % Extract the X and Y positions for phage i at time step k
            x_val = xB_over_time(k, (j-1)*2 + 1) / micron;
            y_val = xB_over_time(k, j*2) / micron;

            fprintf(bacteria_file, '%.6f \t %.6f \t', x_val, y_val);
        end
        fprintf(bacteria_file, '\n');
    end
    fclose(bacteria_file); 
    %fprintf('Trajectories saved to %s\n', filename);

    %% Save Cluster Positions to a .dat File
    clusters_filename = fullfile(outputFolder, 'trajectoriesClusters.dat'); 
    clusters_file = fopen(clusters_filename, 'w'); 

    % Write header for clusters file
    fprintf(clusters_file, 'Time (s)\t');
    for k = 1:num_clusters
        fprintf(clusters_file, 'C%d X (micron) \t C%d Y (micron) \t', k, k);
    end
    fprintf(clusters_file, '\n');

    % Write the bacteria positions at each time step
    for k = 1:N_steps
        fprintf(clusters_file, '%.6f \t', (k-1)*dt); % Time
        for jj = 1:num_clusters
            % Extract the X and Y positions for phage i at time step k
            x_val = xC_over_time(k, (jj-1)*2 + 1) / micron;
            y_val = xC_over_time(k, jj*2) / micron;

            fprintf(clusters_file, '%.6f \t %.6f \t', x_val, y_val);
        end
        fprintf(clusters_file, '\n');
    end
    fclose(clusters_file); 
    %fprintf('Trajectories saved to %s\n', filename);
end
