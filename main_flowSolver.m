clc; clear all; close all;
%% Add paths
addpath('fluid/solverNS/')
addpath('interpolation')
addpath_tools()
%% TODOs:
disp('--------------------------------------------------------------------')
disp('TODO')
disp(' - why is singular all periodic BCs?')
disp('--------------------------------------------------------------------')
%% PARAMETERS: (move to function?)
parameters.case_name  = 'mucus1';
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
minDarcyNum = 1e-3;  %1e-3; % perme->0      : more resistance % adimensional, Darcy number
maxDarcyNum = 1e14; % perme->infty  : no resistance   % adimensional, Darcy number
%parameters.pressureGrad = 0.5*1/minDarcyNum ;       % -20 gradP gives aprox 83.3*1e-6 u_mag
% -----Solver:
parameters.is_parall  = false;
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
%% Check compatibility conditions
if(BC.periodic.topBot && BC.slidingTopBot )
    error('not possible periodicTopBot and slidingTopBot')
end
if(BC.periodic.leftRight && BC.dirichletLeft )
    error('not possible periodic.leftRight and dirichletLeft')
end
if(BC.periodic.topBot==false && BC.periodic.leftRight)
    error('Does not make much sense...')
end
if(BC.dirichletTopBot)
    error('I dont want this BC, too restrictive')
end
%% Mesh
disp('---------------------  Mesh and domain  ----------------------------')
do_mesh_from_image = true;
if(do_mesh_from_image)
    meshName = [parameters.case_name '_mesh'];
    load(['./output/' meshName])
    fprintf('Num nodes: %d\n',size(mesh.X,1))
    
    % Set domain limits:
    domain.x_LL = min(mesh.X(:,1));
    domain.x_LR = max(mesh.X(:,1));
    domain.y_LD = min(mesh.X(:,2));
    domain.y_LU = max(mesh.X(:,2));
    mesh.domain = domain;

    fprintf('Domain:      [%4.2e %4.2e]x[%4.2e %4.2e]\n',domain.x_LL,domain.x_LR,domain.y_LD,domain.y_LU)
 
    % Set permeability:
    perme = mesh.perme; % from image

    min_p = min(perme);
    max_p = max(perme);
    factPerme = (maxDarcyNum-minDarcyNum)/(max_p-min_p);
    perme = minDarcyNum + (perme-min_p)*factPerme;
    model.perme = perme;
    clear perme;
else
    parameters.h_ini = 0.02;
    domain.LLcorner = [-1;-1];
    domain.URcorner = [1 ; 1];
    [mesh] = generateMesh(domain,parameters);

    %model.perme      = @perme_analytical;
    model.perme      = perme_analytical(mesh.X);
end
%% Adimensionalize equation
if(parameters.adimensionalize)
    [domain,mesh,model,parameters]=adimensionalize(domain,mesh,model,parameters);
end
%% Periodicity
if(BC.periodic.leftRight || BC.periodic.topBot)
    checkPeriodicity(mesh);    
    mesh = computePeriodicMaps(mesh);  
end
if(BC.periodic.p.lagrange) 
    model.source     = @(X)(source_zero(X,parameters));
end 
%% Incompressible flow solver
disp('---------------------  Solve flow ----------------------------------')
parameters.case_name = [parameters.case_name '_' BC.label];

mesh = MeshMINI(mesh); 
parameters.model = model;
[solution] = NavierStokesMini(mesh,parameters);

opt_out = solution;
opt_out.name = ['./output/paraview/' parameters.case_name '_sol'];
opt_out.k    = model.perme;
opt_out.m    = solution.node_marks;
exportMeshParaviewSolver(mesh,opt_out);

checkPeriodicSolution(mesh,solution,BC);
%% Build interp
disp('---------------------  Interpolant ---------------------------------')
if(parameters.adimensionalize)
    solution.u_adim = solution.u;
    solution.u = parameters.umag_in_dimensional*solution.u;
end
disp('Setting interpolant...')
[U_interp]=set_interpolant(mesh,solution.u,parameters.adimensionalize,domain);

save(['./output/' parameters.case_name '_uInterp'],"U_interp");
disp('DONE')
%% ALREADY DONE THINGS
% disp(' - adimensionalize equation by characteristic length')
% disp('   - adimensionalize at the beggining and solve stokes')
% disp('   - export for particles with the initial dimensions')
% disp('   -> ADIMENSIONALIZATION DONE')
% disp(' - write adim on report')
% disp('   -> DONE')
% disp(' - program periodic boundary conditions')
% disp('   - mapping of periodic dofs')
% disp('   - solve velocity for unique periodic dofs')
% disp('   - solve pressure with jump on the periodic dofs')
% disp('   - prescribe differential in pressure')