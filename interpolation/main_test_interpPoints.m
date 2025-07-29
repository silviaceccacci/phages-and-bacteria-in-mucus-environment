clear all; close all;
%%
addpath('interpolation/')
addpath('tools/paraview/')
%%
global machine;
machine = 'abelm2';
%%
domain.LLcorner = [-1;-1];
domain.URcorner = [1 ; 1];

parameters.h_ini = 0.1;
[mesh1] = generateMesh(domain,parameters);
parameters.h_ini = 0.2;
[mesh2] = generateMesh(domain,parameters);
parameters.h_ini = 0.05;
[mesh3] = generateMesh(domain,parameters);
dobleDomain.LLcorner = domain.LLcorner*2;
dobleDomain.URcorner = domain.URcorner*2;
parameters.h_ini = 0.025;
[mesh4] = generateMesh(dobleDomain,parameters);

f1 = sin(2*pi*mesh1.X(:,1)).*cos(2*pi*mesh1.X(:,2));
options.exportName = [ 'testInterp_points_1' ];
options.f = f1;
exportMeshParaview(mesh1.X,mesh1.T,options)

% adimensionalize
mesh1_original = mesh1;
domain.Lchar=2; 
domain.origin = domain.LLcorner;  

% domain.Lchar=1;
% domain.origin = 0*domain.LLcorner;  

mesh1.X(:,1) = mesh1.X(:,1)-domain.origin(1);
mesh1.X(:,2) = mesh1.X(:,2)-domain.origin(2);
mesh1.X = mesh1.X/domain.Lchar;
domain.x_LL = min(mesh1.X(:,1));
domain.x_LR = max(mesh1.X(:,1));
domain.y_LD = min(mesh1.X(:,2));
domain.y_LU = max(mesh1.X(:,2));
domain.Lx = domain.x_LR-domain.x_LL;
domain.Ly = domain.y_LU-domain.y_LD;
domain

% domain.Lx = 1;
% domain.Ly = 1;

adim = true;
[F_interp]=set_interpolant(mesh1,f1,adim,domain);
[f2]=interpolatePoints(F_interp,mesh2.X);
[f3]=interpolatePoints(F_interp,mesh3.X);
[f4]=interpolatePoints(F_interp,mesh4.X);
[f5]=interpolatePoints(F_interp,mesh1_original.X);

options.exportName = [ 'testInterp_points_3' ];
options.f = f2;
exportMeshParaview(mesh2.X,mesh2.T,options)

options.exportName = [ 'testInterp_points_2' ];
options.f = f3;
exportMeshParaview(mesh3.X,mesh3.T,options)

options.exportName = [ 'testInterp_points_4' ];
options.f = f4;
exportMeshParaview(mesh4.X,mesh4.T,options)

options.exportName = [ 'testInterp_points_5' ];
options.f = f5;
exportMeshParaview(mesh1_original.X,mesh1_original.T,options)
