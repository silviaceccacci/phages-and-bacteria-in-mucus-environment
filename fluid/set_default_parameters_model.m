function [parameters,model]=set_default_parameters_model()
% ----- PDE model:
parameters.convective      = false;
parameters.adimensionalize = true;
% -----PDE BCs:
[BC_dirichletLeft,BC_lagrangeJump,BC_pgrad,BC_pgrad_noperx]=setDiffernetBCs();
%BC = BC_dirichletLeft; % too artificial, no periodicity on leftRight
%BC = BC_pgrad_noperx;  % periodicity with pressure gradient on equation, no per on press x
%BC = BC_lagrangeJump;  % generalized periodicity using lagrange multi to impose press jump
BC = BC_pgrad;         % periodicity with pressure gradient on equation 
parameters.BC = BC; 
% -----PDE parameters:
model.viscosity         = 1e-2      ; % DISAPEARS IF ADIMENSIONALIZED in V, not in P
parameters.umag_in      = 83.3*1e-6 ; % 
parameters.pressureGrad = -20 ;       % -20 gradP gives aprox 83.3*1e-6 u_mag
parameters.minDarcyNum = 1e-2;  %1e-3; % perme->0      : more resistance % adimensional, Darcy number
parameters.maxDarcyNum = 1e14; % perme->infty  : no resistance   % adimensional, Darcy number
%parameters.pressureGrad = 0.5*1/minDarcyNum ;       % -20 gradP gives aprox 83.3*1e-6 u_mag
% -----Solver:
parameters.is_parall  = true;
if(parameters.is_parall) 
    pool = gcp('nocreate');
    if isempty(pool)
        numWorkers = 4;
        parpool(numWorkers);
    end
end 
parameters.TOL_solver = 1e-3; % for non-linear (Navier-Stokes)
% -----Source term and dirichlet
model.u_dir  = @(X)(bc_const(X,parameters)); 
model.source = @(X)([-parameters.pressureGrad; 0.0]);