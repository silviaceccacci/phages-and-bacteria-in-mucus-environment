function F_val = brinkman_objective(x, case_name, SI_ux_stokes, SI_uy_stokes, TR_stokes, cache_ref)
% BRINKMAN_OBJECTIVE  L2 loss between Brinkman and Stokes velocity fields.
%
%   x(1) = thresh       (solid fraction threshold, e.g. 0.40–0.46)
%   x(2) = log10(Da)    (Darcy number in log scale, e.g. -8 to -6)
%
%   cache_ref is a 1-element cell array holding a struct with fields:
%       .thresh     — threshold used in the last mesh generation
%       .mesh_data  — the loaded mesh .mat struct from that call

thresh = x(1);
Da     = 10^x(2);

% --- Retrieve and update cache -------------------------------------------
cache = cache_ref{1};

if isnan(cache.thresh) || abs(thresh - cache.thresh) > 1e-9
    fprintf('  [mesh]  thresh=%.4f — remeshing...\n', thresh);
    func_main_im2mesh(case_name, thresh);                          % writes mesh to ./output/
    cache.thresh    = thresh;
    cache.mesh_data = load(['./output/' case_name '_mesh.mat']);
    cache_ref{1}    = cache;                                       % write back
end

X_brink_phys = cache.mesh_data.mesh.X;

% --- Map Stokes solution onto Brinkman nodes ------------------------------
u_stokes_at_brink = [SI_ux_stokes(X_brink_phys(:,1), X_brink_phys(:,2)), ...
                     SI_uy_stokes(X_brink_phys(:,1), X_brink_phys(:,2))];

elem_idx = pointLocation(TR_stokes, X_brink_phys(:,1), X_brink_phys(:,2));
in_solid = isnan(elem_idx);
u_stokes_at_brink(in_solid, :) = 0;

% --- Solve Brinkman -------------------------------------------------------
u_brink = func_main_flowSolverBrinkman(case_name, Da, thresh);

% --- L2 loss --------------------------------------------------------------
diff2 = sum((u_brink - u_stokes_at_brink).^2, 2);
norm2 = sum(u_stokes_at_brink.^2, 2);
F_val = sum(diff2) / sum(norm2);

fprintf('  thresh=%.4f | log10(Da)=%6.2f | F=%10.4e\n', thresh, log10(Da), F_val);
end