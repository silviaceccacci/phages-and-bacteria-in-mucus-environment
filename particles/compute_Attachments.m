function [phages, bacteria, attached_phages] = compute_Attachments(phages, bacteria, clusters, lB, d_enc1, d_enc2, d_enc3, last_time_step, outputFolder)
   
     %% Compute bacteria-phage attachment

    % Initialize an array to keep track of attached phages
    num_phages = length(phages);
    num_bacteria = length(bacteria);
    attached_phages = false(num_phages, 1);
    attachment_info = cell(num_bacteria, 1); % Cell array to store attachment details

    % Loop over each phage and each bacterium to check for attachment
    for i = 1:num_phages
        if ~attached_phages(i) % Only check unattached phages
            for j = 1:num_bacteria
                %bact = bacteria(j);
                %distance = norm(phages(i).position - bacteria(j).position);
                distance = pointToSegmentDistance(phages(i).position, bacteria(j).position, bacteria(j).orientation, bacteria(j).length);
                isAttached = (distance < d_enc1);
                %isAttached = check_PhageToRodAttachment(phages(i).position, bacteria(j).position, lB, d_enc1);

                if isAttached
                    phages(i).position = bacteria(j).position;
                    phages(i).velocity = bacteria(j).velocity;
                    phages(i).is_attached = true;
                    bacteria(j).phages_ids(end+1) = phages(i).id;
                    phages(i).attached_bacterium = bacteria(j).id;

                    % Record attachment info
                    if isempty(attachment_info{j})
                        attachment_info{j} = i;
                    else
                        attachment_info{j} = [attachment_info{j}, i];
                    end

                    attached_phages(i) = true;  % Mark the phage as attached
                    break;  % Break the loop to avoid reassigning the phage to another bacterium
                end
            end
        end
    end

    % Only print and save attachment info at the last time step
    if last_time_step
        filename = fullfile(outputFolder, 'attachmentsPhagesBacteria.dat'); 
        fileID = fopen(filename, 'w');

        for j = 1:num_bacteria
            if ~isempty(attachment_info{j})
                fprintf(fileID, '%d', j); % Write bacterium index
                fprintf(fileID, ', %d', attachment_info{j}); % Write attached phages
                fprintf(fileID, '\n');
            end
        end
        fclose(fileID);
        %fprintf('Attachment data saved to %s\n', filename);
    end

    %% Compute bacteriaInCluster-phage attachment

     % Initialize an array to keep track of attached phages
    num_phages = length(phages);
    num_bacteria = length(bacteria);
    attached_phages = false(num_phages, 1);
    attachment_info = cell(num_bacteria, 1); % Cell array to store attachment details

    for i = 1:num_phages
        if ~attached_phages(i) % Only check unattached phages
            for c = 1:length(clusters)
                cluster = clusters(c);
                bacteriaInCluster = cluster.bacteria;

                % Compute distances to all bacteria in the cluster
                distances = arrayfun(@(bact) norm(phages(i).position - bact.position), bacteriaInCluster);

                % Find minimum distance and corresponding bacterium
                [min_dist, idx] = min(distances);

                if min_dist < d_enc2
                    % Attach phage to that bacterium
                    phages(i).position = bacteriaInCluster(idx).position;
                    phages(i).velocity = bacteriaInCluster(idx).velocity;
                    phages(i).is_attached = true;
                    phages(i).attached_bacterium = idx;

                    % Record attachment info
                    if isempty(attachment_info{c})
                        attachment_info{c} = i;
                    else
                        attachment_info{c} = [attachment_info{c}, i];
                    end

                    attached_phages(i) = true;
                    break;  % Stop checking other clusters for this phage
                end
            end
        end
    end

    if last_time_step
        filename = fullfile(outputFolder, 'attachmentsPhagesBacteriaInClusters.dat'); 
        fileID = fopen(filename, 'w');

        for j = 1:num_bacteria
            if ~isempty(attachment_info{j})
                fprintf(fileID, '%d', j); % Write bacterium index
                fprintf(fileID, ', %d', attachment_info{j}); % Write attached phages
                fprintf(fileID, '\n');
            end
        end
        fclose(fileID);
    end

    %% Compute COMCluster-phage attachment

     % Initialize an array to keep track of attached phages
    num_phages = length(phages);
    num_bacteria = length(bacteria);
    attached_phages = false(num_phages, 1);
    attachment_info = cell(length(clusters), 1); % Cell array to store attachment details

    for i = 1:num_phages
        if ~attached_phages(i)  
            for c = 1:length(clusters)
                clusterCOM = clusters(c).position;   % Assuming cluster has a .position field (COM)
                clusterVel = clusters(c).velocity;   % Assuming cluster has a .velocity field

                dist_to_COM = norm(phages(i).position - clusterCOM);

                if dist_to_COM < d_enc3
                    % Attach phage to this cluster's COM
                    phages(i).position = clusterCOM;
                    phages(i).velocity = clusterVel;
                    phages(i).is_attached = true;
                    phages(i).attached_cluster = c;  

                    attached_phages(i) = true;
                    break;  % Stop after first successful attachment
                end
            end
        end
    end

     if last_time_step
        filename = fullfile(outputFolder, 'attachmentsPhagesClustersCOM.dat'); 
        fileID = fopen(filename, 'w');

        for j = 1:length(clusters)
            if ~isempty(attachment_info{j})
                fprintf(fileID, '%d', j); % Write cluster index
                fprintf(fileID, ', %d', attachment_info{j}); % Write attached phages
                fprintf(fileID, '\n');
            end
        end
        fclose(fileID);
     end

end
