function threshold = compute_clustering_threshold(U_interp, x_min, x_max, y_min, y_max, ...
                                                lB, wB, dt, D)
% computeClusteringThreshold
% Computes a reasonable clustering distance for bacteria based on:
%   - bacterium size (length, width)
%   - timestep & diffusion
%   - flow field (u_max)
%   - domain size
%
% Inputs:
%   U_interp : function handle returning [u_x, u_y] at given [x,y]
%   x_min, x_max, y_min, y_max : domain limits
%   lB, wB : bacterium length and width (m)
%   dt : simulation timestep (s)
%   D  : diffusion coefficient (m^2/s)
%
% Output:
%   threshold : clustering distance (m)

% Domain size
Omega_X = x_max - x_min;
Omega_Y = y_max - y_min;

% Grid for sampling velocity
nplot = 100;
hplot = (x_max - x_min) / nplot;
[XX, YY] = meshgrid(x_min:hplot:x_max, y_min:hplot:y_max);

% Evaluate velocity field
UU = zeros(size(XX,1), size(XX,2), 2);
for i = 1:size(XX,1)
    for j = 1:size(XX,2)
        UU(i,j,:) = U_interp([XX(i,j), YY(i,j)]);
    end
end

% Velocity magnitude
UU_mag = sqrt(sum(UU.^2, 3));
u_max = max(UU_mag(:));

%u_max = 83.3e-6; %enforce Da=\infty

% Compute clustering distance
threshold = compute_clustering_distance(lB, wB, dt, D, u_max, Omega_X, Omega_Y);

disp(['Max velocity: ', num2str(u_max), ' m/s'])
disp(['Clustering distance threshold: ', num2str(threshold*1e6), ' μm'])

end
