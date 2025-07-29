clear all; close all;
%%
addpath('interpolation/')
addpath('tools/paraview/')

disp('to execute, move to main folder')
%%
global machine;
machine = 'abelm2';
%%
domain.LLcorner = [-1;-1];
domain.URcorner = [1 ; 1];

parameters.h_ini = 0.1;
[mesh1] = generateMesh(domain,parameters);

f1 = sin(2*pi*mesh1.X(:,1)).*cos(2*pi*mesh1.X(:,2));

parameters.h_ini = 0.2;
[mesh2] = generateMesh(domain,parameters);

[f2]=interpolateMeshes(mesh1,mesh2,f1);

options.exportName = [ 'testInterp_1' ];
options.f = f1;
exportMeshParaview(mesh1.X,mesh1.T,options)

options.exportName = [ 'testInterp_2' ];
options.f = f2;
exportMeshParaview(mesh2.X,mesh2.T,options)


domain.LLcorner = domain.LLcorner*2;
domain.URcorner = domain.URcorner*2;
parameters.h_ini = 0.1;
[mesh3] = generateMesh(domain,parameters);
[f3]=interpolateMeshes(mesh1,mesh3,f1);

options.exportName = [ 'testInterp_3' ];
options.f = f3;
exportMeshParaview(mesh3.X,mesh3.T,options)
