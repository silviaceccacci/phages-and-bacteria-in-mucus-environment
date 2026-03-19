function [phages, bacteria, attached_phages, n_phages_attached] = compute_attachments_3(phages, bacteria, clusters, max_num_bacteria, attached_phages, ...
    d_enc1, d_enc2, d_enc3, last_time_step, outputFolder)

n_phages_attached = 0;
attachment_info = cell(max_num_bacteria, 1);

phages_attached = find(attached_phages);
phages_not_attached = find(~attached_phages);

if isempty(phages_not_attached)
    return;
end

if ~isempty(bacteria)
    positions_bacteria = cat(1, bacteria.position);
else
    positions_bacteria = [];
end

if ~isempty(clusters)
    positions_COM = cat(1, clusters.position);
else
    positions_COM = [];
end

positions_phages = cat(1, phages(phages_not_attached).position);

n_phages = length(phages_not_attached);
n_clusters = length(clusters);

if ~isempty(positions_bacteria)
    D_bacteria = pdist2(positions_bacteria, positions_phages);
    [min_D_bacteria, id_bacteria] = min(D_bacteria, [], 1);
else
    min_D_bacteria = inf(1, n_phages);
    id_bacteria = zeros(1, n_phages);
end

if ~isempty(clusters)
    D_bacteria_in_cluster = inf(n_clusters, n_phages);
    idbac_bacteria_in_cluster = zeros(n_clusters, n_phages);
    for local_i = 1:n_clusters
        bacterias_cluster = clusters(local_i).bacteria;
        if isempty(bacterias_cluster)
            continue;
        end
        positions_bacteria_in_cluster = cat(1, bacterias_cluster.position);
        [dist, iddd] = min(pdist2(positions_bacteria_in_cluster, positions_phages));
        D_bacteria_in_cluster(local_i,:) = dist;
        idbac_bacteria_in_cluster(local_i,:) = iddd;
    end
    [min_D_bacteria_in_cluster, id_bacteria_in_cluster] = min(D_bacteria_in_cluster, [], 1);
else
    min_D_bacteria_in_cluster = inf(1, n_phages);
    id_bacteria_in_cluster = zeros(1, n_phages);
end

if ~isempty(positions_COM)
    D_COM = pdist2(positions_COM, positions_phages);
    [min_D_COM, id_COM] = min(D_COM, [], 1);
else
    min_D_COM = inf(1, n_phages);
    id_COM = zeros(1, n_phages);
end

A = [min_D_bacteria./d_enc1; min_D_bacteria_in_cluster./d_enc2; min_D_COM./d_enc3];

[min_D, id_particle] = min(A, [], 1);
idx_phage_to_particle = find(min_D < 1);
id_particle = id_particle(idx_phage_to_particle);

%Attach phages to single bacteria
if ~isempty(bacteria)
    phages_to_bacteria = idx_phage_to_particle(id_particle == 1);
    for k = 1:length(phages_to_bacteria)
        local_i = phages_to_bacteria(k);
        global_i = phages_not_attached(local_i);
        j = id_bacteria(local_i);
        if j < 1 || j > numel(bacteria)
            continue; % skip invalid index
        end
        n_phages_attached = n_phages_attached + 1;
        phages(global_i).position = bacteria(j).position;
        phages(global_i).velocity = bacteria(j).velocity;
        phages(global_i).is_attached = true;
        phages(global_i).attached_bacterium = bacteria(j).id;

        if ~ismember(phages(global_i).id, bacteria(j).phages_ids)
            bacteria(j).phages_ids(end+1) = phages(global_i).id;
        end

        if isempty(attachment_info{j})
            attachment_info{j} = global_i;
        else
            attachment_info{j} = [attachment_info{j}, global_i];
        end
        attached_phages(global_i) = true;
    end
end

%Attach phages to bacteria inside clusters
if ~isempty(clusters)
    phages_to_bacteria_cluster = idx_phage_to_particle(id_particle == 2);
    for k = 1:length(phages_to_bacteria_cluster)
        local_i = phages_to_bacteria_cluster(k);
        global_i = phages_not_attached(local_i);
        c = id_bacteria_in_cluster(local_i);
        if c < 1 || c > numel(clusters)
            continue;
        end
        cluster = clusters(c);
        bacteriaInCluster = cluster.bacteria;
        if isempty(bacteriaInCluster)
            continue;
        end
        j = idbac_bacteria_in_cluster(c, local_i);
        if j < 1 || j > numel(bacteriaInCluster)
            continue;
        end
        n_phages_attached = n_phages_attached + 1;
        phages(global_i).position = bacteriaInCluster(j).position;
        phages(global_i).velocity = bacteriaInCluster(j).velocity;
        phages(global_i).is_attached = true;
        phages(global_i).attached_bacterium = bacteriaInCluster(j).id;
        attached_phages(global_i) = true;
    end
end

%Attach phages to cluster COM 
if ~isempty(clusters)
    phages_to_COM = idx_phage_to_particle(id_particle == 3);
    for k = 1:length(phages_to_COM)
        local_i = phages_to_COM(k);
        global_i = phages_not_attached(local_i);
        j = id_COM(local_i);
        if j < 1 || j > numel(clusters)
            continue;
        end
        phages(global_i).position = clusters(j).position;
        phages(global_i).velocity = clusters(j).velocity;
        phages(global_i).is_attached = true;
        attached_phages(global_i) = true;
    end
end

%Update already attached phages
for k = 1:length(phages_attached)
    local_i = phages_attached(k);
    id_of_bacteria_with_phage_i = phages(local_i).attached_bacterium;
    id_of_cluster_with_phage_i = phages(local_i).attached_cluster;

    if id_of_bacteria_with_phage_i ~= -1 && ~isempty(bacteria)
        [idx_bacteria, idx_cluster, in_cluster] = ...
            index_of_bacteria_with_id(clusters, bacteria, id_of_bacteria_with_phage_i);
        if in_cluster == 0
            phages(local_i).position = bacteria(idx_bacteria).position;
            phages(local_i).velocity = bacteria(idx_bacteria).velocity;
        elseif in_cluster == 1
            phages(local_i).position = clusters(idx_cluster).bacteria(idx_bacteria).position;
            phages(local_i).velocity = clusters(idx_cluster).bacteria(idx_bacteria).velocity;
        end
    end

    if id_of_cluster_with_phage_i ~= -1 && ~isempty(clusters)
        [idx_COM] = find_closest_COM_to_phage(clusters, phages(local_i));
        phages(local_i).position = clusters(idx_COM).position;
        phages(local_i).velocity = clusters(idx_COM).velocity;
    end
end

end
