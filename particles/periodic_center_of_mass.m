function x_com = periodic_center_of_mass(x, L)
    % x: vector of particle positions
    % L: domain length (e.g., 2 for [-1,1] or [0,2])
    theta = 2*pi*x/L;
    C = mean(cos(theta));
    S = mean(sin(theta));
    theta_com = atan2(S, C);     % mean angle
    x_com = mod(theta_com * L/(2*pi), L);
end