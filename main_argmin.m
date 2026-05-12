clc; clear all; close all;

%% Add paths
addpath('./input')
addpath('./image2mesh/src')
addpath('fluid/')
addpath('fluid/solverNS/')
addpath('interpolation')
addpath('argmin')
addpath_tools()

%% Parameters
case_name = 'PSD_topography';
fixed_stokes_threshold = 0.7; % Baseline threshold for solid Stokes mesh

nThresh = 4;
nDa = 2;
thresh_frac_vec = linspace(0.40, 0.46, nThresh);
Da_vec = logspace(-8, -6, nDa);

F = nan(nThresh, nDa);
rho_solid = nan(nThresh, 1);

%% =========================================================================
%  STEP 1 — Stokes Reference (Fixed)
% =========================================================================
fprintf('=== STEP 1: Generating Stokes Reference Solution ===\n');

% 1. Create Stokes solid mesh (saves to ./output/...)
func_Fluid_Solid_mesh(case_name, fixed_stokes_threshold);

% 2. Solve Stokes
[u_stokes, X_stokes_phys, T_stokes] = func_main_flowSolverStokes(case_name);

% 3. Build scattered interpolant using PHYSICAL coordinates
SI_ux_stokes = scatteredInterpolant(X_stokes_phys(:,1), X_stokes_phys(:,2), u_stokes(:,1), 'linear', 'none');
SI_uy_stokes = scatteredInterpolant(X_stokes_phys(:,1), X_stokes_phys(:,2), u_stokes(:,2), 'linear', 'none');
TR_stokes = triangulation(T_stokes, X_stokes_phys(:,1), X_stokes_phys(:,2));

%% =========================================================================
%  STEP 2 — Parameter Sweep
% =========================================================================
fprintf('\n=== STEP 2: Sweep %d thresholds x %d Da = %d solves ===\n', nThresh, nDa, nThresh*nDa);

for i = 1:nThresh
    thresh = thresh_frac_vec(i);
    fprintf('\n--- Threshold %d/%d (frac=%.3f) ---\n', i, nThresh, thresh);

    % 4. Create Brinkman porous mesh (saves to ./output/...)
    rho_solid(i) = func_main_im2mesh(case_name, thresh);

    % 5. Load the Brinkman mesh just created to get physical coordinates
    temp = load(['./output/' case_name '_mesh.mat']);
    X_brink_phys = temp.mesh.X;

    current_mesh = temp.mesh; 
    rho_current = rho_solid(i);

    % 6. Map Stokes onto Brinkman nodes (Done once per threshold)
    u_stokes_at_brink = [SI_ux_stokes(X_brink_phys(:,1), X_brink_phys(:,2)), ...
                         SI_uy_stokes(X_brink_phys(:,1), X_brink_phys(:,2))];

    % Zero out nodes inside solid obstacles
    elem_idx = pointLocation(TR_stokes, X_brink_phys(:,1), X_brink_phys(:,2));
    in_solid = isnan(elem_idx);
    u_stokes_at_brink(in_solid, :) = 0;

    % 7. Inner loop for Darcy parameter
    for j = 1:nDa
        Da = Da_vec(j);

        % Solve Brinkman
        [u_brink] = func_main_flowSolverBrinkman(case_name, Da, thresh);
        save_filename = sprintf('./output/%s_Brinkman_Da%.2e_thresh%.2f.mat', case_name, Da, thresh);
        save(save_filename, 'current_mesh', 'u_brink', 'Da', 'rho_current', 'thresh');

        % Error minimization (L2 Difference)
        diff2 = sum((u_brink - u_stokes_at_brink).^2, 2);
        F(i,j) = mean(diff2);
        
        fprintf('  Da=%8.2e | F=%10.4e\n', Da, F(i,j));
    end
end

%% =========================================================================
%  STEP 3 — Find Argmin and Plot
% =========================================================================
F_clean = F;
F_clean(isnan(F)) = inf;

[F_min, idx] = min(F_clean(:));
[i_opt, j_opt] = ind2sub([nThresh, nDa], idx);
thresh_opt = thresh_frac_vec(i_opt);
Da_opt = Da_vec(j_opt);
rho_opt = rho_solid(i_opt);

fprintf('\n=====================================================\n');
fprintf('  ARGMIN RESULT\n');
fprintf('  threshold* = %.3f  (solid fraction rho* = %.3f)\n', thresh_opt, rho_opt);
fprintf('  Da* = %.4e\n', Da_opt);
fprintf('  F* = %.4e\n', F_min);
fprintf('=====================================================\n\n');

% Plotting (Surface)
[DA_grid, TH_grid] = meshgrid(log10(Da_vec), thresh_frac_vec);
F_plot = log10(F + 1e-30);

figure('Position', [100 100 900 600]);
surf(DA_grid, TH_grid, F_plot, 'EdgeAlpha', 0.2); colormap(turbo); shading interp; colorbar;
hold on;
plot3(log10(Da_opt), thresh_opt, log10(F_min + 1e-30), 'rp', 'MarkerSize', 24, 'MarkerFaceColor', 'red');
xlabel('log_{10}(Da_{min})'); ylabel('threshold fraction'); zlabel('log_{10}(F)');
title('Argmin: F(threshold, Da) = ||u_{Brinkman} - u^*_{Stokes}||^2');
grid on; view([-35, 30]);

save('./output/argmin_results_modular.mat', 'F', 'thresh_frac_vec', 'Da_vec', 'rho_solid', 'thresh_opt', 'Da_opt', 'F_min');

%% =========================================================================
%  STEP 4 — Visual Comparison (Normalized Contour Plots)
% =========================================================================
fprintf('\n=== STEP 4: Generating Normalized Contour Plots ===\n');

% 1. Load the best-fit Brinkman solution we saved to disk
best_filename = sprintf('./output/%s_Brinkman_Da%.2e_thresh%.2f.mat', case_name, Da_opt, thresh_opt);
fprintf('Loading best-fit data from: %s\n', best_filename);
best_data = load(best_filename);
X_brink_opt = best_data.current_mesh.X;
u_brink_opt = best_data.u_brink;

% 2. Calculate Velocity Magnitudes
mag_stokes = sqrt(u_stokes(:,1).^2 + u_stokes(:,2).^2);
mag_brink  = sqrt(u_brink_opt(:,1).^2 + u_brink_opt(:,2).^2);

% 3. Normalize the magnitudes (0 to 1 scale) so patterns can be compared
mag_stokes_norm = mag_stokes / max(mag_stokes);
mag_brink_norm  = mag_brink / max(mag_brink);

% 4. Create a unified, dense regular grid for contour plotting
num_pts = 300; % Resolution of the contour plot
xq = linspace(min(X_stokes_phys(:,1)), max(X_stokes_phys(:,1)), num_pts);
yq = linspace(min(X_stokes_phys(:,2)), max(X_stokes_phys(:,2)), num_pts);
[Xq, Yq] = meshgrid(xq, yq);

% 5. Interpolate both unstructured meshes onto the regular grid
SI_stokes = scatteredInterpolant(X_stokes_phys(:,1), X_stokes_phys(:,2), mag_stokes_norm, 'linear', 'none');
Vq_stokes = SI_stokes(Xq, Yq);

SI_brink = scatteredInterpolant(X_brink_opt(:,1), X_brink_opt(:,2), mag_brink_norm, 'linear', 'none');
Vq_brink = SI_brink(Xq, Yq);

% 6. Calculate the Absolute Difference
Vq_diff = abs(Vq_stokes - Vq_brink);

% Find valid pixels (ignore NaN values outside the mesh boundaries)
valid_mask = ~isnan(Vq_stokes) & ~isnan(Vq_brink);

% 1. Correlation Coefficient (R) -> 1.0 is a perfect structural match
R_matrix = corrcoef(Vq_stokes(valid_mask), Vq_brink(valid_mask));
R_value  = R_matrix(1,2); 

% 2. Mean Absolute Error (MAE) -> Average pixel difference (0.0 is perfect)
MAE_norm = mean(Vq_diff(valid_mask));

fprintf('--- Contour Similarity Metrics ---\n');
fprintf('Correlation (R): %.4f\n', R_value);
fprintf('Mean Abs Error : %.4f\n\n', MAE_norm);

% 7. Plotting (Only the Difference)
figure('Position', [200, 200, 600, 500], 'Name', 'Absolute Difference Contour');
contourf(Xq, Yq, Vq_diff, 30, 'LineColor', 'none');
colormap(hot); % 'hot' colormap makes error heatmaps look great
colorbar; 
clim([0, max(max(Vq_diff(:)), 0.01)]);

title(sprintf('Absolute Difference |Stokes - Brinkman|\n(Da* = %.1e, \\rho* = %.3f) | Pattern Match (R) = %.3f', Da_opt, rho_opt, R_value), 'FontSize', 13);
xlabel('x', 'FontSize', 12); 
ylabel('y', 'FontSize', 12);
axis equal tight;