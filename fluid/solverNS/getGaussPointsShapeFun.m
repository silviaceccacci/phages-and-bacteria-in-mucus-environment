function [N,Nxi,Neta,NP,NxiP,NetaP,pospg,wpg,numel,nen,nenP,numnp,ngaus,ndofn,nunk,nunkP,ngeom]=...
    getGaussPointsShapeFun(elem,IEN,X,IENP,XP)

[numel,nen] = size(IEN); 
[numelP,nenP] = size(IENP);
% Total number of nodes
if elem == 1 && nen == 4 % MINI
    numnp = size(X,1) + numel;
else
    numnp = size(X,1);     
end
% Degrees of freedom
nunk  = 2*numnp; ndofn = 2*nen;
nunkP = size(XP,1);

% Number of Gauss points used in the numerical quadrature (ngaus)
% and number of nodes to interpolate the geometry (ngeom)
if elem == 0
    if nen == 4
        ngaus = 4; ngeom = 4;
    elseif nen == 9
        ngaus = 9; ngeom = 9;
    end
else
    if nen == 3
        ngaus = 3; ngeom = 3;
    elseif nen == 6 
        ngaus = 4; ngeom = 6;
    elseif nen == 4 
        ngaus = 7; ngeom = 3;
    end
end

% Quadrature and shape functions on the Gauss points
[pospg,wpg] = Quadrature(elem,ngaus);
[N,Nxi,Neta] = ShapeFunc(elem,nen,pospg);
[NP,NxiP,NetaP] = ShapeFunc(elem,nenP,pospg);

