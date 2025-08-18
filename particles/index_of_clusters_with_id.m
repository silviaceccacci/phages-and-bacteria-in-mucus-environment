function [idx_COM] = index_of_clusters_with_id(clusters, id_of_interest)
for i = 1 : length(clusters)
    id = clusters(i).id;
    if id == id_of_interest
        idx_COM = i;
        str = sprintf('index of cluster con fago in com = %d', idx_COM);
        disp(str)
    else
        disp('id cluster mai uguale a id of interest')
   end
end
end