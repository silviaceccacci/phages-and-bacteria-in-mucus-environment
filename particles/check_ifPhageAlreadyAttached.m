function phages = check_ifPhageAlreadyAttached(phages, bacteria)
    num_phages = length(phages);
    num_bacteria = length(bacteria);
    for i = 1:num_phages
        if isfield(phages(i), "is_attached") && phages(i).is_attached
            b_idx = phages(i).attached_bacterium;
            if b_idx > 0 && b_idx <= num_bacteria
                phages(i).position = bacteria(b_idx).position;
                phages(i).velocity = bacteria(b_idx).velocity;
            else
                warning('Phage %d refers to invalid bacterium index %d', i, b_idx);
            end
        end
    end
end