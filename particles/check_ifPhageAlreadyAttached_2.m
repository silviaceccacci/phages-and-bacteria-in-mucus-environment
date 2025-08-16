function [phages, attached_phages] = check_ifPhageAlreadyAttached_2(phages, bacteria, attached_phages)
    for i = 1:length(phages)
        disp('ENTRO NELLA FUNCTION check_ifPhageAlreadyAttached_2')
        %if isfield(phages(i), "is_attached") && phages(i).is_attached
        if phages(i).is_attached
            disp('SONO UN FAGO ATTACHED')
            b_idx = phages(i).attached_bacterium;
            if b_idx > 0 && b_idx <= length(bacteria)
                str = sprintf('ID del batterio a cui sono attached = %d', b_idx);
                disp(str);
                %phages(i).relative_position_wrt_bact = phages(i).position - bacteria(b_idx).position;
                %rel_pos = phages(i).relative_position_wrt_bact;
                %phages(i).position = bacteria(b_idx).position + rel_pos;
                phages(i).position = bacteria(b_idx).position;
                phages(i).velocity = bacteria(b_idx).velocity;
                attached_phages(i) = true;
            else
                warning('Phage %d refers to invalid bacterium index %d', i, b_idx);
            end
        end
    end
end