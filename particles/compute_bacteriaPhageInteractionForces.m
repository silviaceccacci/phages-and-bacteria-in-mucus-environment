function [bacteriumInteractionForces, phageInteractionForces] = compute_bacteriaPhageInteractionForces_2(bacteria, phages, epsilon,phages_radius)
    
    % Initialize bacterium and phage force accumulators for bacterium-phage interaction
    num_bacteria = length(bacteria);
    num_phages = length(phages);
    bacteriumInteractionForces = zeros(num_bacteria, 2);
    phageInteractionForces = zeros(num_phages, 2);

    % Calculate bacterium-phage interaction forces
    phages_position = cat(1, phages(:).position);
    if(nargin==3) 
        phages_radius = cat(1, phages(:).radius);
    end
    for j = 1:num_bacteria
        vij = phages_position - bacteria(j).position;
        r = sqrt(sum(vij.^2, 2));
        sigma = 2^(-1/6) * (bacteria(j).width + phages_radius); % Distance for bacterium-phage interaction
        r_cutoff = 2^(1/6) * sigma; % Cutoff distance
        idx = find(r < r_cutoff);
        sr6=(sigma ./ r).^6;
        %F_magnitude = 24 * epsilon * (-2 * (sigma ./ r).^12 + (sigma ./ r).^6);
        F_magnitude = 24 * epsilon * (-2 * sr6.^2 + sr6);
        F_interaction = zeros(num_phages, 2);
        F_interaction(idx, :) = F_magnitude(idx) .* (vij(idx, :)./r(idx)); 
        bacteriumInteractionForces(j, :) = bacteriumInteractionForces(j, :) + sum(F_interaction, 1);
        phageInteractionForces = phageInteractionForces - F_interaction;

%         for i = 1:num_phages
%             % Calculate vector from phage to bacterium
%             rij = abs(bacteria(j).radius - phages(i).radius);
%             vij = norm(rij); % Distance magnitude
%             sigma = 2^(-1/6) * (bacteria(j).width + phages(i).radius); % Distance for bacterium-phage interaction
%             r_cutoff = 2^(1/6) * sigma; % Cutoff distance
% 
%             if vij < r_cutoff
%                 % Calculate the WCA potential (not needed but useful to remember)
%                 %V = 4 * epsilon * ((sigma / r)^12 - (sigma / r)^6 + epsilon);
% 
%                 % Calculate the magnitude of the force
%                 F_magnitude = 24 * epsilon * (-2 * (sigma / vij)^12 + (sigma / vij)^6);
% 
%                 % Calculate the force vector
%                 F_interaction = F_magnitude * (rij / vij); % Force is in the direction of rij
% 
%                 %F_interaction = [0, 0];
%             else
%                 % Beyond the cutoff distance, potential and force are zero
%                 %V = 0;
%                 F_interaction = [0, 0];
%             end
% 
%             % Accumulate forces on bacterium and phage
%             bacteriumInteractionForces(j, :) = bacteriumInteractionForces(j, :) + F_interaction;
%             phageInteractionForces(i, :) = phageInteractionForces(i, :) - F_interaction; % Equal and opposite force on phage
%         end
    end
end
