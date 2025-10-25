function [phages, bacteria, attached_phages, n_phages_attached] = compute_attachments(phages, bacteria, clusters, max_num_bacteria, attached_phages, d_enc1, d_enc2, d_enc3, last_time_step, outputFolder)

n_phages_attached = 0;
attachment_info = cell(max_num_bacteria, 1);
%     str = sprintf('num bacteria in fun1 = %d', num_bacteria);
%     disp(str);

for i = 1:length(phages) % Loop over each phage and each bacterium to check for attachment
    flag = 0;
    if ~attached_phages(i) % Only check unattached phages
        %check attachment to single bacterium
        for j = 1:length(bacteria)
            distance = norm(phages(i).position - bacteria(j).position);
            isAttached = (distance < d_enc1);

            if isAttached
                n_phages_attached = n_phages_attached + 1;
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
                %disp('I should be attached in un batterio singolo')
                flag = 1;
                break;  % Break the loop to avoid reassigning the phage to another bacterium
            end
        end

        %check attachment on bacteria in clusters
        if flag == 0
            for c = 1:length(clusters)
                cluster = clusters(c);
                bacteriaInCluster = cluster.bacteria;
    
                % Compute distances to all bacteria in the cluster
                distances = arrayfun(@(bact) norm(phages(i).position - bact.position), bacteriaInCluster);
    
                % Find minimum distance and corresponding bacterium
                [min_dist, idx] = min(distances);
    
                if min_dist < d_enc2
                    n_phages_attached = n_phages_attached + 1;
                    %phages(i).relative_position_wrt_bact = phages(i).position - bacteriaInCluster(idx).position;
                    %rel_pos = phages(i).relative_position_wrt_bact;
                    %phages(i).position = bacteriaInCluster(idx).position + rel_pos;
                    phages(i).position = bacteriaInCluster(idx).position;
                    %disp('Questa e la prima posizione del fago')
                    %phages(i).position
    
                    phages(i).velocity = bacteriaInCluster(idx).velocity;
                    phages(i).is_attached = true;
                    phages(i).attached_bacterium = bacteriaInCluster(idx).id;
                    %str = sprintf('id batterio nel cluster da funzione min_dist = %d', bacteriaInCluster(idx).id);
                    %disp(str);
    
                    if isempty(attachment_info{c})
                        attachment_info{c} = i;
                    else
                        attachment_info{c} = [attachment_info{c}, i];
                    end
    
                    attached_phages(i) = true;
                    %disp('I should be attached in un batterio nel cluster')
                    flag = 1;
                    break;
                end
    
            end
        end

        %check for attachment to cluster COM
        if flag == 0
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
                    %phages(i).attached_cluster = clusters(c).id;
                    %str = sprintf('id cluster COM = %d', clusters(c).id);
                    %disp(str);
    
                    if isempty(attachment_info{c})
                        attachment_info{c} = i;
                    else
                        attachment_info{c} = [attachment_info{c}, i];
                    end
    
                    attached_phages(i) = true;
                    flag = 1;
                    break;
                end
            end
        end


    else
        id_of_bacteria_with_phage_i = phages(i).attached_bacterium;
        %str = sprintf('id batterio con fago = %d', id_of_bacteria_with_phage_i);
        %disp(str);
        id_of_cluster_with_phage_i = phages(i).attached_cluster;
        %str = sprintf('id cluster con fago in com = %d', id_of_cluster_with_phage_i);
        %disp(str);
        if id_of_bacteria_with_phage_i ~= -1
            %disp('sono nella id del batterio')
            [idx_bacteria, idx_cluster, in_cluster] = index_of_bacteria_with_id(clusters, bacteria, id_of_bacteria_with_phage_i);
            if in_cluster == 0
                %disp('il batterio sta da solo' )
                phages(i).position = bacteria(idx_bacteria).position;
                phages(i).velocity = bacteria(idx_bacteria).velocity;
            end
            if in_cluster == 1
                %disp('il batterio sta nel cluster')
                phages(i).position = clusters(idx_cluster).bacteria(idx_bacteria).position;
                phages(i).velocity = clusters(idx_cluster).bacteria(idx_bacteria).velocity;
                %disp('Questa e la posizione del fago dopo primo attachment')
                %phages(i).position
            end

        end
        if id_of_cluster_with_phage_i ~= -1
            %disp('sono nel COM del cluster')
            %[idx_COM] = index_of_clusters_with_id(clusters, id_of_cluster_with_phage_i);
            [idx_COM] = find_closest_COM_to_phage(clusters, phages(i));
            phages(i).position = clusters(idx_COM).position;
            %str = sprintf('position of COM = %d', clusters(idx_COM).position);
            %disp(str);
            phages(i).velocity = clusters(idx_COM).velocity;

        end

    end
end

% if last_time_step % Only print and save attachment info at the last time step
%     filename = fullfile(outputFolder, 'attachmentsPhagesBacteria.dat');
%     fileID = fopen(filename, 'w');
% 
%     %         str = sprintf('ciao ');
%     %         disp(str);
% 
%     for j = 1:length(bacteria)
%         if ~isempty(attachment_info{j})
%             fprintf(fileID, '%d', bacteria(j).id);
%             %                 str = sprintf('bacteria ID = %d', bacteria(j).id);
%             %                 disp(str);
%             fprintf(fileID, ', %d', attachment_info{j}); % Write attached phages
%             fprintf(fileID, '\n');
%         end
%     end
%     fclose(fileID);
%     %fprintf('Attachment data saved to %s\n', filename);
% end


end
