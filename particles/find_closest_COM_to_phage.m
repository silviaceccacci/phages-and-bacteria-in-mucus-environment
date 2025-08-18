function idx_cluster_COM_clostest_to_phage_of_interest = find_closest_COM_to_phage(clusters, phage_of_interest)
distance_min = Inf;
for i = 1:length(clusters)
    distance = norm(clusters(i).position - phage_of_interest.position);
    if distance < distance_min
        distance_min = distance;
        idx_cluster_COM_clostest_to_phage_of_interest = i;
    end
end
end
