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
parameters.case_name    = 'mucus1';
parameters.umag_in      = 83.3*1e-6; % from some reference
parameters.pressureGrad = -20  ; % -20 gradP gives aprox 83.3*1e-6 u_mag for some perme
parameters.minDarcyNum  = 1e-5 ; % perme->0      : more resistance % adimensional, Darcy number
parameters.maxDarcyNum  = 1e14 ; % perme->infty  : no resistance   % adimensional, Darcy number
%% Mesh
disp('---------------------  Mesh and domain  ----------------------------')
do_mesh_from_image  = true;
[mesh,model,domain] = generate_mesh_domain(do_mesh_from_image,model,parameters);
%% Flow solver
umag_mucus    = parameters.umag_in;
if(parameters.convective) 
    % nonlinear equation requires Newton to rescale gradPress for umean
    maxIterVel    = 10;
    isConvMeanVel = false;
    fact_pressGrad= 1;
    iter_vel      = 0;
    tol_u         = 1e-2;
    while(isConvMeanVel==false)
        iter_vel = iter_vel+1;
        parameters.pressureGrad = parameters.pressureGrad*fact_pressGrad;
    
        [solution,domain]=flow_solver(domain,mesh,model,parameters);
        
        umean = computeMeanVel(mesh,domain,solution);
        fact_pressGrad = umag_mucus/umean;
        
        isConvMeanVel = abs(fact_pressGrad-1)<tol_u || iter_vel>maxIterVel;
        %isConvMeanVel = abs(umag_mucus-umean)<tol_u*umag_mucus && iter_vel<maxIterVel;
        
        disp('----------------------------------------------------------------------------------- ')
        fprintf(' Target mucus vel  |  Mean flow vel  | Factor pGrad |  Press Grad |  isConverged \n')
        fprintf('      %8.2e     |    %8.2e     |    %6.2f    |    %6.2f   |    %s\n', ...
        umag_mucus, umean, fact_pressGrad, parameters.pressureGrad, string(isConvMeanVel));
        disp('----------------------------------------------------------------------------------- ')
    end
else    
    % Linear equation: equation linear on gradPress, rescale u & p directly
    % For the stokes problem, due to the liniarity of the problem, the
    % desired mean velocity is attained by rescaling the solution
    [solution,domain]=flow_solver(domain,mesh,model,parameters);
    umean            = computeMeanVel(mesh,domain,solution);
    fact_rescale     = umag_mucus/umean;
    solution.u = solution.u*fact_rescale;
    solution.p = solution.p*fact_rescale;
    umean            = computeMeanVel(mesh,domain,solution);
    fprintf('Mean velocity: %8.2e\n',umean)
end
%% Build interp
disp('---------------------  Interpolant ---------------------------------')
disp('Setting interpolant...')
[U_interp]=set_interpolant(mesh,solution.u,false,domain);

save(['./output/' parameters.case_name '_uInterp'],"U_interp");
