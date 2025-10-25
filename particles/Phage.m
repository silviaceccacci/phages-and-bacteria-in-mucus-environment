classdef Phage
    properties
        id
        radius
        mass
        friction_coefficient
        noise_amplitude
        position
        velocity
        domain
        x_min
        y_min
        x_max
        y_max
        is_attached
        attached_bacterium
        attached_cluster
        relative_position_wrt_bact
        relative_position_wrt_com
    end
    
    methods
        function p = Phage(rP, rhoP, mu_water, kB, T, domain, x_min, y_min, x_max, y_max)
            p.radius = rP;
            p.mass = rhoP * 4/3 * pi * rP^3;
            p.friction_coefficient = 6 * pi * mu_water * rP;
            p.noise_amplitude = sqrt(2 * p.friction_coefficient * kB * T); 
            p.position = rand(1, 2) .* domain; % Random initial position
            p.velocity = zeros(1, 2); % Initial velocity zero
            p.domain = domain;
            p.x_min = x_min;
            p.y_min = y_min;
            p.x_max = x_max;
            p.y_max = y_max;
            p.is_attached = false;
            p.attached_bacterium = -1; %-1 if unattached
            p.attached_cluster = -1; %-1 if unattached
        end
        
        function [p, noise_term, force] = computeFluidForce(p, fluid_velocity)
            % Compute the fluid drag force on the phage
            noise_term = p.noise_amplitude * randn(1, 2);
            force = -p.friction_coefficient * (p.velocity - fluid_velocity) + noise_term;
        end
        
        function p = updateVelocity(p, force, dt)
            % Update velocity using Langevin equation
            p.velocity = p.velocity + (force / p.mass) * dt;
        end
        
        function p = updatePosition(p, dt, domain, x_min, y_min, x_max, y_max)
            % Update position based on velocity and ensure boundary conditions
            p.position = p.position + p.velocity * dt;
            
%             % Reflect at boundaries if needed
%             if p.position(1) < 0 || p.position(1) > domain(1)
%                 p.velocity(1) = -p.velocity(1); % Reverse x-velocity
%             end
%             if p.position(2) < 0 || p.position(2) > domain(2)
%                 p.velocity(2) = -p.velocity(2); % Reverse y-velocity
%             end

            % Periodic boundary conditions
            if p.position(1) < x_min
                p.position(1) = p.position(1) + domain(1); 
            elseif p.position(1) > x_max
                p.position(1) = p.position(1) - domain(1); 
            end
            if p.position(2) < y_min
                p.position(2) = p.position(2) + domain(2); 
            elseif p.position(2) > y_max
                p.position(2) = p.position(2) - domain(2); 
            end
        end

        function p = attachTo(p, bacteriumIndex)
            p.is_attached = true;
            p.attached_bacterium = bacteriumIndex;
        end

        function p = detach(p)
            p.is_attached = false;
            p.attached_bacterium = -1;
        end
    end
end
