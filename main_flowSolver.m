clc; clear all; close all;
%% Add paths
addpath('fluid/')
addpath('fluid/solverNS/')
addpath('interpolation')
addpath_tools()
%% TODOs:
disp('--------------------------------------------------------------------')
disp('TODO')
disp(' - iterate on flow to set proper pressure for mean flow')
disp('--------------------------------------------------------------------')
%% PARAMETERS: (move to function?)
[parameters,model]=set_default_parameters_model();
parameters.case_name  = 'mucus1';
parameters.umag_in      = 83.3*1e-6 ; % from some reference
parameters.pressureGrad = -20 ; % -20 gradP gives aprox 83.3*1e-6 u_mag for some perme
parameters.minDarcyNum = 1e-2;  %1e-3; % perme->0      : more resistance % adimensional, Darcy number
parameters.maxDarcyNum = 1e14;  %perme->infty  : no resistance   % adimensional, Darcy number
%% Mesh
disp('---------------------  Mesh and domain  ----------------------------')
do_mesh_from_image = true;
[mesh,model,domain]=generate_mesh_domain(do_mesh_from_image,model,parameters);
%% Flow solverrameters)
[solution,domain]=flow_solver(domain,mesh,model,parameters);
%% Build interp
disp('---------------------  Interpolant ---------------------------------')
disp('Setting interpolant...')
[U_interp]=set_interpolant(mesh,solution.u,parameters.adimensionalize,domain);

save(['./output/' parameters.case_name '_uInterp'],"U_interp");
