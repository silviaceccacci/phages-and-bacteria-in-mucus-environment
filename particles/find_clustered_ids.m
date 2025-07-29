function clustered_ids = find_clustered_ids(bacteria, clustered_bacteria)
    clustered_ids = [];
    for b = clustered_bacteria
        for i = 1:length(bacteria)
            if norm(bacteria(i).position - b.position) < 1e-10
                clustered_ids(end+1) = i;
                break;
            end
        end
    end
end