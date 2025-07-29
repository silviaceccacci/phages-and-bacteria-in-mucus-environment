function save_trajectories(uP_over_time, uB_over_time, dt, N_steps, num_phages, num_bacteria, micron, outputFolder)
    % Function to save phage and bacteria velocities to .dat files

     % Ensure the output folder exists
    if ~exist(outputFolder, 'dir')
        mkdir(outputFolder);
    end

    % Save Phage Positions to a .dat File
    %phage_filename = 'phages_trajectories.dat';
    phage_filename = fullfile(outputFolder, 'velocitiesPhages.dat'); 
    phage_file = fopen(phage_filename, 'w'); 

    % Write header for phage file
    fprintf(phage_file, 'Time (s) \t');
    for i = 1:num_phages
        fprintf(phage_file, 'P%d ux (m/s) \t P%d uy (m/s) \t', i, i);
        %fprintf(phage_file, 'Phage %d X (micron) \t Phage %d Y (micron) \t', i, i);
    end
    fprintf(phage_file, '\n');

    % Write the phage positions at each time step
    for k = 1:N_steps
        fprintf(phage_file, '%.6f \t', (k-1)*dt); % Time
        for i = 1:num_phages
            % Extract the X and Y positions for phage i at time step k
            x_val = uP_over_time(k, (i-1)*2 + 1) / micron;
            y_val = uP_over_time(k, i*2) / micron;

            % Debugging: print values to console to verify correctness
            %fprintf('Time: %.6f, Phage %d X: %.6f, Y: %.6f\n', (k-1)*dt, i, x_val, y_val);
            
            fprintf(phage_file, '%.6f \t %.6f \t', x_val, y_val);
        end
        fprintf(phage_file, '\n');
    end
    fclose(phage_file); 

    
    % Save Bacteria Positions to a .dat File
    %bacteria_filename = 'bacteria_trajectories.dat';
    bacteria_filename = fullfile(outputFolder, 'velocitiesBacteria.dat'); % Updated path
    bacteria_file = fopen(bacteria_filename, 'w'); 

    % Write header for bacteria file
    fprintf(bacteria_file, 'Time (s)\t');
    for j = 1:num_bacteria
        fprintf(bacteria_file, 'B%d ux (m/s) \t B%d uy (m/s) \t', j, j);
        %fprintf(bacteria_file, 'Bacteria %d X (micron) \t Bacteria %d Y (micron) \t', j, j);
    end
    fprintf(bacteria_file, '\n');

    % Write the bacteria positions at each time step
    for k = 1:N_steps
        fprintf(bacteria_file, '%.6f \t', (k-1)*dt); % Time
        for j = 1:num_bacteria
            % Extract the X and Y positions for phage i at time step k
            x_val = uB_over_time(k, (j-1)*2 + 1) / micron;
            y_val = uB_over_time(k, j*2) / micron;

            fprintf(bacteria_file, '%.6f \t %.6f \t', x_val, y_val);
        end
        fprintf(bacteria_file, '\n');
    end
    fclose(bacteria_file); 
    %fprintf('Velocities saved to %s\n', filename);
end
