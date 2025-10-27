function d_attach = compute_clustering_distance(L, W, dt, D, u_max, Omega_X, Omega_Y)
% computeClusteringDistance
% ---------------------------------------------------------------
% Computes a reasonable clustering distance for rod-shaped bacteria
% in a given domain, accounting for geometry, motion, and domain size.
%
% Inputs:
%   L       : bacterium length (m)
%   W       : bacterium width (m)
%   dt      : simulation timestep (s)
%   D       : diffusion coefficient (m^2/s)
%   u_max   : max flow velocity (m/s)
%   Omega_X : domain size x (m)
%   Omega_Y : domain size y (m)
%
% Output:
%   d_attach : clustering distance (m)
% ---------------------------------------------------------------

% 1) Geometric contribution (edges + small fraction of length)
fraction_of_length = 0.15;  % tune between 0.1-0.25
d_geom = W + fraction_of_length*L;

% 2) Motion contribution (diffusion + flow)
d_motion = sqrt(2*D*dt) + u_max*dt;

% 3) Combine
d_attach_raw = d_geom + d_motion;

% 4) Cap by fraction of domain to avoid excessive clustering
max_fraction_of_domain = 0.15;  % 10-20% of smallest domain dimension
d_max = max_fraction_of_domain * min(Omega_X, Omega_Y);

d_attach = min(d_attach_raw, d_max);
end
