 classdef Cluster

    properties
        id
        position
        velocity
        bacteria
        friction_coefficient
        domain
        x_min
        y_min
        x_max
        y_max
        bacteria_ids
    end

    properties (Dependent)
        size %computed as number of bacteria
    end

    methods
        function c = Cluster(position, velocity, bacteria, domain, x_min, y_min, x_max, y_max)
            if nargin > 0
                c.position = rand(1, 2);
                c.velocity = rand(1, 2);
                c.bacteria = bacteria;
                c.friction_coefficient = 0;
                c.domain = domain;
                c.x_min = x_min;
                c.y_min = y_min;
                c.x_max = x_max;
                c.y_max = y_max;
            end
        end

        function N = get.size(c)
            N = length(c.bacteria);
        end

        function c = update_friction_coefficient(c, mu, wB)
            N = c.size;
            effective_radius = wB * N^(1/2); %N^(1/3) in 3D
            c.friction_coefficient = 6 * pi * mu * effective_radius;
        end

        function c = update_center_of_mass_and_group_velocity(c)
            num_bacteria = length(c.bacteria);
            masses = zeros(num_bacteria, 1);
            positions = zeros(num_bacteria, 2);
            velocities = zeros(num_bacteria, 2);

            for i = 1:num_bacteria
                masses(i) = c.bacteria(i).mass;
                positions(i,:) = c.bacteria(i).position;
                velocities(i,:) = c.bacteria(i).velocity;
            end

            total_mass = sum(masses);
            center_of_mass = sum(positions.*masses, 1)/total_mass;
            c.position = center_of_mass;
            group_velocity = sum(velocities.*masses, 1)/total_mass;
            c.velocity = group_velocity;

            for i = 1:num_bacteria
                %c.bacteria(i).relative_velocity = c.bacteria(i).velocity - c.velocity;
                c.bacteria(i).relative_velocity = c.velocity;
            end

            for i = 1:num_bacteria
                c.bacteria(i).relative_position_wrt_com = c.bacteria(i).position - c.position;
            end
        end

        function c = evolveLangevin(c, fluid_velocity, dt, kB, T, domain, x_min, y_min, x_max, y_max)
            noise = sqrt(2 * c.friction_coefficient * kB * T / dt) * randn(1, 2);
            drag_force = -c.friction_coefficient * (c.velocity - fluid_velocity);
            total_force = drag_force + noise;

            num_bacteria = length(c.bacteria);
            masses = zeros(num_bacteria, 1);
            for i = 1:num_bacteria
                masses(i) = c.bacteria(i).mass;
            end
            total_mass = sum(masses);

            c.velocity = c.velocity + (total_force ./ total_mass) * dt;
            c.position = c.position + c.velocity * dt;

            for i = 1:c.size
                rel_pos = c.bacteria(i).relative_position_wrt_com;
                %rel_pos = c.bacteria(i).position - c.position; % relative to old COM
                c.bacteria(i).position = c.position + rel_pos; % shifted with new COM
                c.bacteria(i).velocity = c.velocity; % move with cluster

                if c.bacteria(i).position(1) < x_min
                    c.bacteria(i).position(1) = c.bacteria(i).position(1) + domain(1);
                elseif c.bacteria(i).position(1) > x_max
                    c.bacteria(i).position(1) = c.bacteria(i).position(1) - domain(1);
                end
                if c.bacteria(i).position(2) < y_min
                    c.bacteria(i).position(2) = c.bacteria(i).position(2) + domain(2);
                elseif c.bacteria(i).position(2) > y_max
                    c.bacteria(i).position(2) = c.bacteria(i).position(2) - domain(2);
                end
            end

            if c.position(1) < x_min
                c.position(1) = c.position(1) + domain(1);
            elseif c.position(1) > x_max
                c.position(1) = c.position(1) - domain(1);
            end
            if c.position(2) < y_min
                c.position(2) = c.position(2) + domain(2);
            elseif c.position(2) > y_max
                c.position(2) = c.position(2) - domain(2);
            end
        end

%             for i = 1:c.size
%                 rel_pos = c.bacteria(i).relative_position_wrt_com;
%                 %rel_pos = c.bacteria(i).position - c.position; % relative to old COM
%                 c.bacteria(i).position = c.position + rel_pos; % shifted with new COM
%                 c.bacteria(i).velocity = c.velocity; % move with cluster
%             end
% 
%             % Check if COM or any bacterium is out of bounds
%             shift_x = 0;
%             shift_y = 0;
% 
%             % Check X direction
%             if c.position(1) < x_min || any(arrayfun(@(b) b.position(1) < x_min, c.bacteria))
%                 shift_x = domain(1); % move right
%             elseif c.position(1) > x_max || any(arrayfun(@(b) b.position(1) > x_max, c.bacteria))
%                 shift_x = -domain(1); % move left
%             end
% 
%             % Check Y direction
%             if c.position(2) < y_min || any(arrayfun(@(b) b.position(2) < y_min, c.bacteria))
%                 shift_y = domain(2); % move up
%             elseif c.position(2) > y_max || any(arrayfun(@(b) b.position(2) > y_max, c.bacteria))
%                 shift_y = -domain(2); % move down
%             end
% 
%             % Apply shift to entire cluster if needed
%             if shift_x ~= 0 || shift_y ~= 0
%                 c.position = c.position + [shift_x, shift_y];
%                 for i = 1:c.size
%                     c.bacteria(i).position = c.bacteria(i).position + [shift_x, shift_y];
%                 end
%             end
%        end
    end

    methods (Static)
        function [clusters, bacteriaNotInCluster] = form_clusters(bacteriaList, threshold, domain, x_min, y_min, x_max, y_max)

            num_bacteria = length(bacteriaList);
            visited = false(1, num_bacteria);
            clusters = [];
            cluster_id_counter = 1;

            for i = 1:num_bacteria
                if visited(i)
                    continue;
                end

                group = i;
                visited(i) = true;
                queue = i;

                while ~isempty(queue)
                    current = queue(1);
                    queue(1) = [];

                    for j = 1:num_bacteria
                        if ~visited(j) %if not visited
                            dist = norm(bacteriaList(current).position - bacteriaList(j).position);
                            if dist < threshold
                                visited(j) = true;
                                group(end+1) = j;
                                queue(end+1) = j;
                            end
                        end
                    end
                end

                if numel(group) > 1

                bacteriaInCluster = bacteriaList(group);
                pos = mean(cat(1, bacteriaInCluster.position), 1);
                vel = mean(cat(1, bacteriaInCluster.velocity), 1);

                clusters = [clusters, Cluster(pos, vel, bacteriaInCluster,domain, x_min, y_min, x_max, y_max)];
                cluster_id_counter = cluster_id_counter + 1;
%                 str = sprintf('cluster_id_counter = %d', cluster_id_counter);
%                 disp(str);
                else
                    visited(i) = false;
                end
            end
            bacteriaNotInCluster = bacteriaList(~visited);
        end
    end
end
