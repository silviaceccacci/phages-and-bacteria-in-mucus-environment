function [mesh] = generateMesh(domain,parameters)
   
doPlotMeshGen = false;

domain.x_LL = domain.LLcorner(1);
domain.x_LR = domain.URcorner(1);
domain.y_LD = domain.LLcorner(2);
domain.y_LU = domain.URcorner(2);

Lx = domain.x_LR-domain.x_LL;
Ly = domain.y_LU-domain.y_LD;
ne1dx = ceil(Lx/parameters.h_ini);
ne1dy = ceil(Ly/parameters.h_ini);

[X,T] = generateRectangleMesh(ne1dx,ne1dy);
% New domain [x_LL,x_LR]x[y_LL,y_LR]
for i=1:size(X,1)
    X(i,1)=((-domain.x_LL+domain.x_LR)/2)*X(i,1)+((-domain.x_LL+domain.x_LR)/2)+domain.x_LL;
    X(i,2)=((-domain.y_LD+domain.y_LU)/2)*X(i,2)+((-domain.y_LD+domain.y_LU)/2)+domain.y_LD;
end

% n1d=100;
% xh = linspace(x_LL,x_LR,n1d);
% yh = linspace(y_LD,y_LU,n1d);
% [xx,yy] = meshgrid(xh,yh);
% X = [X; xx(:) yy(:)];
% 
% T = delaunay(X);

mesh.domain = domain;
mesh.X = X;
mesh.T = T;

% if(doPlotMeshGen)
%     options.exportName = [parameters.case_name '_mesh'];
%     exportMeshParaview(X,T,options);
% end
    