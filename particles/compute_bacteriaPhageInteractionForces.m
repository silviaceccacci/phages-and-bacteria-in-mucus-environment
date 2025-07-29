function [bacteriumInteractionForces, phageInteractionForces] = compute_bacteriaPhageInteractionForces(bacteria, phages, epsilon)
    
    % Initialize bacterium and phage force accumulators for bacterium-phage interaction
    num_bacteria = length(bacteria);
    num_phages = length(phages);
    bacteriumInteractionForces = zeros(num_bacteria, 2);
    phageInteractionForces = zeros(num_phages, 2);

    % Calculate bacterium-phage interaction forces
    for j = 1:num_bacteria
        for i = 1:num_phages
            % Calculate vector from phage to bacterium
            rij = abs(bacteria(j).radius - phages(i).radius);
            r = norm(rij); % Distance magnitude
            sigma = 2^(-1/6) * (bacteria(j).width + phages(i).radius); % Distance for bacterium-phage interaction
            r_cutoff = 2^(1/6) * sigma; % Cutoff distance

            if r < r_cutoff
                % Calculate the WCA potential (not needed but useful to remember)
                %V = 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6 + epsilon);

                % Calculate the magnitude of the force
                F_magnitude = 24 * epsilon * (-2 * (sigma / r)^12 + (sigma / r)^6);

                % Calculate the force vector
                F_interaction = F_magnitude * (rij / r); % Force is in the direction of rij
            else
                % Beyond the cutoff distance, potential and force are zero
                %V = 0;
                F_interaction = [0, 0];
            end

            % Accumulate forces on bacterium and phage
            bacteriumInteractionForces(j, :) = bacteriumInteractionForces(j, :) + F_interaction;
            phageInteractionForces(i, :) = phageInteractionForces(i, :) - F_interaction; % Equal and opposite force on phage
        end
    end
end
