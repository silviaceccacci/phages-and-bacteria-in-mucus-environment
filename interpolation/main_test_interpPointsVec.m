clear all; close all;
%%
addpath('interpolation/')
addpath('image2mesh/src/paraview/')
%%
global machine;
machine = 'abelm2';
%%
domain.LLcorner = [-1;-1];
domain.URcorner = [1 ; 1];

parameters.h_ini = 0.1;
[mesh1] = generateMesh(domain,parameters);

fv1 = sin(2*pi*mesh1.X(:,1));
fv2 = cos(2*pi*mesh1.X(:,2));
f1 = [fv1,fv2];

[F_interp]=set_interpolant(mesh1,f1);
%F_interp would be given as an ouput of the flow solution
% And now we want to evaluate it on out fagos

parameters.h_ini = 0.2;
[mesh2] = generateMesh(domain,parameters);
[f2]=interpolatePoints(F_interp,mesh2.X);

parameters.h_ini = 0.05;
[mesh3] = generateMesh(domain,parameters);
[f3]=interpolatePoints(F_interp,mesh3.X);

options.exportName = [ 'testInterp_pVEC_1' ];
options.f = f1;
exportMeshParaview(mesh1.X,mesh1.T,options)

options.exportName = [ 'testInterp_pVEC_2' ];
options.f = f2;
exportMeshParaview(mesh2.X,mesh2.T,options)

options.exportName = [ 'testInterp_pVEC_3' ];
options.f = f3;
exportMeshParaview(mesh3.X,mesh3.T,options)