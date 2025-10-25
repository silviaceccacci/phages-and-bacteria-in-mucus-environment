classdef Bacterium
    properties
        id
        length
        width
        mass
        friction_coefficient
        noise_amplitude
        run_velocity
        tumble_prob
        position
        velocity
        phi
        flagellar_force
        fluid_force
        domain
        x_min
        y_min
        x_max
        y_max
        relative_velocity
        relative_position_wrt_com
        phages_ids
        radius
    end
    
    methods
        function b = Bacterium(lB, wB, rhoB, mu_water, kB, T, dt, run_velocity, tumble_frequency, domain, x_min, y_min, x_max, y_max)
            b.length = lB;
            b.width = wB;
            b.mass = rhoB * pi * lB * wB^2;
            b.friction_coefficient = 6 * pi * mu_water * wB;
            b.noise_amplitude = sqrt(2 * b.friction_coefficient * kB * T); %divide by dt?
            b.run_velocity = run_velocity;
            b.tumble_prob = tumble_frequency * dt; % Probability of a tumble in each time step
            b.position = rand(1, 2) .* domain; % Random initial position
            b.velocity = zeros(1, 2); % Initial velocity zero
            b.phi = 0;
            b.flagellar_force = zeros(1, 2);
            b.fluid_force = zeros(1, 2);
            b.domain = domain;
            b.x_min = x_min;
            b.y_min = y_min;
            b.x_max = x_max;
            b.y_max = y_max;
            b.radius = 1e-6;
        end
        
        function b = computePropulsionForce(b)
            % Decide whether to tumble or run
            if rand < b.tumble_prob
                % Tumble phase: stop the motion
                b.velocity = [0, 0];
                b.flagellar_force = [0, 0];
            else
                % Run phase: compute random angles for flagellar force
                b.phi = 2 * pi * rand; 
                propulsion_dir = [cos(b.phi),sin(b.phi)];
                b.flagellar_force = b.friction_coefficient * b.run_velocity * propulsion_dir;
            end
        end

        function b = computeFluidForce(b, fluid_velocity)
            % Compute the fluid drag force on the bacterium
            noise_term = b.noise_amplitude * randn(1, 2);
            b.fluid_force = -b.friction_coefficient * (b.velocity - fluid_velocity) + noise_term;
        end

        function total_force = computeTotalBacteriumForce(b)
            % Compute the total force acting on the bacterium
            total_force = b.flagellar_force + b.fluid_force;
        end
        
        function b = updateVelocity(b, total_force, dt)
            % Update velocity based on total force
            b.velocity = b.velocity + (total_force / b.mass) * dt;
        end
        
        function b = updatePosition(b, dt, domain, x_min, y_min, x_max, y_max)
            b.position = b.position + b.velocity * dt;
            
%             % Reflective boundary conditions 
%             if b.position(1) < x_min || b.position(1) > x_max
%                 b.velocity(1) = -b.velocity(1); 
%             end
%             if b.position(2) < y_min || b.position(2) > y_max
%                 b.velocity(2) = -b.velocity(2); 
%             end

            % Periodic boundary conditions
            if b.position(1) < x_min
                b.position(1) = b.position(1) + domain(1); 
            elseif b.position(1) > x_max
                b.position(1) = b.position(1) - domain(1); 
            end
            if b.position(2) < y_min
                b.position(2) = b.position(2) + domain(2); 
            elseif b.position(2) > y_max
                b.position(2) = b.position(2) - domain(2); 
            end
        end
    end
end
