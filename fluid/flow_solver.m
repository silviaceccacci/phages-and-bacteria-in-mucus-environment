function [solution,domain]=flow_solver(domain,mesh,model,parameters)
%% Check compatibility conditions
BC = parameters.BC;
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

%% Dimensionalize solution
if(parameters.adimensionalize)
    solution.u_adim = solution.u;
    solution.u = parameters.umag_in_dimensional*solution.u;
end