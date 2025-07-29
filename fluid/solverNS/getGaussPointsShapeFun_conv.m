function [N,Nxi,Neta,pospg,pespg,numel,nen,numnp,ngaus,nnodes,ndofn,nunk]=...
    getGaussPointsShapeFun_conv(elem,IEN,X)

% Total number of elements and element node's number
[numel,nen] = size(IEN); 
% Total number of nodes
if elem == 1 && nen == 4
    numnp = size(X,1) + numel;
else
    numnp = size(X,1);     
end
% Degrees of freedom
ndofn = 2*nen; nunk  = 2*numnp; 

% Number of Gauss points used in the numerical quadrature (ngaus)
% and number of nodes to interpolate the geometry (ngeom)
if elem == 0
    if nen == 4
        ngaus = 4; nnodes = 4;
    elseif nen == 9
        ngaus = 9; nnodes = 9;
    end
else
   if nen == 3
        ngaus = 7; nnodes = 3;
   elseif nen == 6 
        ngaus = 7; nnodes = 6;
   elseif nen == 4 
       ngaus = 12; nnodes = 3;
   end
end

% Quadrature and shape functions on the Gauss points
[pospg,pespg] = Quadrature(elem,ngaus);
[N,Nxi,Neta] = ShapeFunc(elem,nen,pospg);
