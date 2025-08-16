function [idx_bacteria,idx_cluster, in_cluster] = index_of_bacteria_with_id(clusters, bacteria, id_of_interest)
idx_bacteria = -10;
idx_cluster = -10;
for i = 1 : length(bacteria)
    id = bacteria(i).id;
    if id == id_of_interest
        idx_bacteria = i;
        in_cluster = 0;
    end
end

if idx_bacteria == -10
    for i = 1 : length(clusters)
        for j = 1 : length(clusters(i).bacteria)
            id = clusters(i).bacteria(j).id;
            if id == id_of_interest
                idx_bacteria = j;
                idx_cluster = i;
                in_cluster = 1;
            end
        end
    end
end

end