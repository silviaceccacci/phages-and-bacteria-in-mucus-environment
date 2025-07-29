function [bacteriumInteractionForces, phageInteractionForces] = compute_bacteriaPhageInteractionForces2(bacteria, phages, epsilon)

    num_bacteria = length(bacteria);
    num_phages = length(phages);
    bacteriumInteractionForces = zeros(num_bacteria, 2);
    phageInteractionForces = zeros(num_phages, 2);

    for j = 1:num_bacteria
        center = bacteria(j).position;              % [x, y]
        %orientation = bacteria(j).orientation;      % unit vector [ux, uy]
        half_length = 0.5 * bacteria(j).length;
        p1 = center - half_length; 
        p2 = center + half_length; 
        %p1 = center - half_length * orientation;    % one end of rod
        %p2 = center + half_length * orientation;    % other end

        for i = 1:num_phages
            q = phages(i).position;                 % phage position (point)

            % Vector from phage to closest point on rod segment
            [r_vec, r] = pointToSegmentVector(q, p1, p2);

            % Effective contact distance
            sigma = 2^(-1/6) * (bacteria(j).width + phages(i).radius);
            r_cutoff = 2^(1/6) * sigma;

            if r < r_cutoff
                % Compute magnitude of WCA force
                F_mag = 24 * epsilon * (-2 * (sigma / r)^12 + (sigma / r)^6);
                F_interaction = F_mag * (r_vec / r);  % direction: phage to rod
            else
                F_interaction = [0, 0];
            end

            % Apply equal and opposite forces
            bacteriumInteractionForces(j, :) = bacteriumInteractionForces(j, :) + F_interaction;
            phageInteractionForces(i, :) = phageInteractionForces(i, :) - F_interaction;
        end
    end
end

function [r_vec, r_mag] = pointToSegmentVector(q, p1, p2)
    v = p2 - p1;
    w = q - p1;
    t = dot(w, v) / dot(v, v);
    t = max(0, min(1, t));  % clamp to [0,1]
    closest = p1 + t * v;
    r_vec = q - closest;
    r_mag = norm(r_vec);
end
