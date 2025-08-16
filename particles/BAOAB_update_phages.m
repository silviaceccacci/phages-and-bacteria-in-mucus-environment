function p = BAOAB_update(p, fluid_velocity, dt, kB, T)
    % BAOAB Langevin integrator for underdamped Langevin dynamics
    % fluid_velocity = local fluid velocity at phage position
    
    m = p.mass;
    gamma = p.friction_coefficient; % friction coefficient
    F = -gamma * (p.velocity - fluid_velocity); % deterministic drag
    
    % === B-step: Half force kick ===
    p.velocity = p.velocity + (0.5 * dt / m) * F;
    
    % === A-step: Half position drift ===
    p.position = p.position + 0.5 * dt * p.velocity;
    
    % === O-step: Ornstein–Uhlenbeck thermalization ===
    exp_fac = exp(-gamma * dt / m);
    noise_std = sqrt((kB * T / m) * (1 - exp_fac^2));
    p.velocity = exp_fac * p.velocity + noise_std * randn(1, 2);
    
    % === A-step: Second half position drift ===
    p.position = p.position + 0.5 * dt * p.velocity;
    
    % === B-step: Second half force kick ===
    % Recompute force at new position
    F_new = -gamma * (p.velocity - fluid_velocity);
    p.velocity = p.velocity + (0.5 * dt / m) * F_new;
    
    % === Periodic BCs ===
    if p.position(1) < p.x_min
        p.position(1) = p.position(1) + p.domain(1);
    elseif p.position(1) > p.x_max
        p.position(1) = p.position(1) - p.domain(1);
    end
    if p.position(2) < p.y_min
        p.position(2) = p.position(2) + p.domain(2);
    elseif p.position(2) > p.y_max
        p.position(2) = p.position(2) - p.domain(2);
    end
end
