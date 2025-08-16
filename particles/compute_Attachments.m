function [phages, bacteria, attached_phages] = compute_attachments(phages, bacteria, clusters, max_num_bacteria, attached_phages, d_enc1, d_enc2, d_enc3, last_time_step, outputFolder)

% %% 0. Always update position of already attached phages
% for i = 1:length(phages)
%     if attached_phages(i)
%         % Attached to a bacterium (case 1 or case 2)
%         if isfield(phages(i), 'attached_bacterium') && ~isempty(phages(i).attached_bacterium)
%             disp('------------------------------------------ ENTROOOOOOO 1111111 ------------------------------------------');
%             %str = sprintf('num bacteria in fun1 = %d', num_bacteria);
%             %disp(str);
%             % First, check if it's a regular bacterium (case 1)
%             idx_bact = find([bacteria.id] == phages(i).attached_bacterium, 1);
%             if ~isempty(idx_bact)
%                 phages(i).position = bacteria(idx_bact).position;
%                 phages(i).velocity = bacteria(idx_bact).velocity;
%                 continue
%             end
% 
%             % Otherwise, it might be a bacterium inside a cluster (case 2)
%             found = false;
%             for c = 1:length(clusters)
%                 disp('------------------------------------------ ENTROOOOOOO 2222222 ------------------------------------------');
%                 bacteriaInCluster = clusters(c).bacteria;
%                 idx_bact_cluster = find([bacteriaInCluster.id] == phages(i).attached_bacterium, 1);
%                 if ~isempty(idx_bact_cluster)
%                     phages(i).position = bacteriaInCluster(idx_bact_cluster).position;
%                     phages(i).velocity = bacteriaInCluster(idx_bact_cluster).velocity;
%                     found = true;
%                     break
%                 end
%             end
%             if found, continue; end
%         end
% 
%         % Attached to a cluster COM (case 3)
%         if isfield(phages(i), 'attached_cluster') && ~isempty(phages(i).attached_cluster)
%             disp('------------------------------------------ ENTROOOOOOO 33333333 ------------------------------------------');
%             c = phages(i).attached_cluster;
%             phages(i).position = clusters(c).position;
%             phages(i).velocity = clusters(c).velocity;
%         end
%     end

%% 1. Compute bacteria-phage attachment
attachment_info = cell(max_num_bacteria, 1);
%     str = sprintf('num bacteria in fun1 = %d', num_bacteria);
%     disp(str);

for i = 1:length(phages) % Loop over each phage and each bacterium to check for attachment
    if ~attached_phages(i) % Only check unattached phages
        for j = 1:length(bacteria)
            distance = norm(phages(i).position - bacteria(j).position);
            isAttached = (distance < d_enc1);

            if isAttached
                %phages(i).relative_position_wrt_bact = phages(i).position - bacteria(j).position; %Delta = x1' - x1
                %rel_pos = phages(i).relative_position_wrt_bact;
                %phages(i).position = bacteria(j).position + rel_pos; %x1' = x1 + Delta
                phages(i).position = bacteria(j).position;

                phages(i).velocity = bacteria(j).velocity;
                phages(i).is_attached = true;

                if ~ismember(phages(i).id, bacteria(j).phages_ids)
                    bacteria(j).phages_ids(end+1) = phages(i).id;
                end
                %bacteria(j).phages_ids(end+1) = phages(i).id;
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

    else
        for j = 1:length(bacteria)
            distance = norm(phages(i).position - bacteria(j).position);
            isAttached = (distance < d_enc1);

            if isAttached
                disp('sono attached e mi muovo con il batterio')
                %phages(i).relative_position_wrt_bact = phages(i).position - bacteria(j).position; %Delta = x1' - x1
                %rel_pos = phages(i).relative_position_wrt_bact;
                %phages(i).position = bacteria(j).position + rel_pos; %x1' = x1 + Delta
                phages(i).position = bacteria(j).position;

                phages(i).velocity = bacteria(j).velocity;
                %phages(i).is_attached = true;

                if ~ismember(phages(i).id, bacteria(j).phages_ids)
                    bacteria(j).phages_ids(end+1) = phages(i).id;
                end
                %bacteria(j).phages_ids(end+1) = phages(i).id;
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

if last_time_step % Only print and save attachment info at the last time step
    filename = fullfile(outputFolder, 'attachmentsPhagesBacteria.dat');
    fileID = fopen(filename, 'w');

    %         str = sprintf('ciao ');
    %         disp(str);

    for j = 1:length(bacteria)
        if ~isempty(attachment_info{j})
            fprintf(fileID, '%d', bacteria(j).id);
            %                 str = sprintf('bacteria ID = %d', bacteria(j).id);
            %                 disp(str);
            fprintf(fileID, ', %d', attachment_info{j}); % Write attached phages
            fprintf(fileID, '\n');
        end
    end
    fclose(fileID);
    %fprintf('Attachment data saved to %s\n', filename);
end

%% 2. Compute bacteriaInCluster-phage attachment
attachment_info = cell(max_num_bacteria, 1);

for i = 1:length(phages)
    if ~attached_phages(i)
        for c = 1:length(clusters)
            cluster = clusters(c);
            bacteriaInCluster = cluster.bacteria;

            % Compute distances to all bacteria in the cluster
            distances = arrayfun(@(bact) norm(phages(i).position - bact.position), bacteriaInCluster);

            % Find minimum distance and corresponding bacterium
            [min_dist, idx] = min(distances);

            if min_dist < d_enc2
                %phages(i).relative_position_wrt_bact = phages(i).position - bacteriaInCluster(idx).position;
                %rel_pos = phages(i).relative_position_wrt_bact;
                %phages(i).position = bacteriaInCluster(idx).position + rel_pos;
                phages(i).position = bacteriaInCluster(idx).position;

                phages(i).velocity = bacteriaInCluster(idx).velocity;
                phages(i).is_attached = true;
                phages(i).attached_bacterium = idx;

                if isempty(attachment_info{c})
                    attachment_info{c} = i;
                else
                    attachment_info{c} = [attachment_info{c}, i];
                end

                attached_phages(i) = true;
                disp('I should be attached')
                break;
            end

        end
    else
        disp('sono nel cluster e mi muovo con il batterio')
        for c = 1:length(clusters)
            cluster = clusters(c);
            bacteriaInCluster = cluster.bacteria;

            % Compute distances to all bacteria in the cluster
            distances = arrayfun(@(bact) norm(phages(i).position - bact.position), bacteriaInCluster);

            % Find minimum distance and corresponding bacterium
            [min_dist, idx] = min(distances);

            if min_dist < d_enc2
                %phages(i).relative_position_wrt_bact = phages(i).position - bacteriaInCluster(idx).position;
                %rel_pos = phages(i).relative_position_wrt_bact;
                %phages(i).position = bacteriaInCluster(idx).position + rel_pos;
                phages(i).position = bacteriaInCluster(idx).position;

                phages(i).velocity = bacteriaInCluster(idx).velocity;
                %phages(i).is_attached = true;
                phages(i).attached_bacterium = idx;

                if isempty(attachment_info{c})
                    attachment_info{c} = i;
                else
                    attachment_info{c} = [attachment_info{c}, i];
                end

                attached_phages(i) = true;
                disp('I should be attached')
                break;
            end

        end
    end
end

if last_time_step
    filename = fullfile(outputFolder, 'attachmentsPhagesBacteriaInClusters.dat');
    fileID = fopen(filename, 'w');

    for j = 1:length(bacteria)
        if ~isempty(attachment_info{j})
            fprintf(fileID, '%d', bacteria(j).id);
            %                 str = sprintf('bacteria ID = %d', bacteria(j).id);
            %                 disp(str);
            fprintf(fileID, ', %d', attachment_info{j});
            fprintf(fileID, '\n');
        end
    end
    fclose(fileID);
end

%% 3. Compute COMCluster-phage attachment
attachment_info = cell(length(clusters), 1);

for i = 1:length(phages)
    if ~attached_phages(i)
        for c = 1:length(clusters)
            clusterCOM = clusters(c).position;
            %phages(i).relative_position_wrt_com = phages(i).position - clusterCOM;
            %rel_pos = phages(i).relative_position_wrt_com;
            clusterVel = clusters(c).velocity;

            dist_to_COM = norm(phages(i).position - clusterCOM);

            if dist_to_COM < d_enc3
                %phages(i).position = clusterCOM + rel_pos;
                phages(i).position = clusterCOM;

                phages(i).velocity = clusterVel;
                phages(i).is_attached = true;
                phages(i).attached_cluster = c;

                if isempty(attachment_info{c})
                    attachment_info{c} = i;
                else
                    attachment_info{c} = [attachment_info{c}, i];
                end

                attached_phages(i) = true;
                break;
            end
        end

    else
        disp('')
        for c = 1:length(clusters)
            clusterCOM = clusters(c).position;
            %phages(i).relative_position_wrt_com = phages(i).position - clusterCOM;
            %rel_pos = phages(i).relative_position_wrt_com;
            clusterVel = clusters(c).velocity;

            dist_to_COM = norm(phages(i).position - clusterCOM);

            if dist_to_COM < d_enc3
                %phages(i).position = clusterCOM + rel_pos;
                phages(i).position = clusterCOM;

                phages(i).velocity = clusterVel;
                phages(i).is_attached = true;
                phages(i).attached_cluster = c;

                if isempty(attachment_info{c})
                    attachment_info{c} = i;
                else
                    attachment_info{c} = [attachment_info{c}, i];
                end

                attached_phages(i) = true;
                break;
            end
        end
    end
end

if last_time_step
    filename = fullfile(outputFolder, 'attachmentsPhagesClustersCOM.dat');
    fileID = fopen(filename, 'w');

    for j = 1:length(clusters)
        if ~isempty(attachment_info{j})
            fprintf(fileID, '%d', clusters(j).id);
            %                 str = sprintf('cluster ID = %d', clusters(j).id);
            %                 disp(str);
            fprintf(fileID, ', %d', attachment_info{j});
            fprintf(fileID, '\n');
        end
    end
    fclose(fileID);
end

end
