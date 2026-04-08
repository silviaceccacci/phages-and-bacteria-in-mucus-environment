function [Ke,Ge,fe] = EleMat(Xe,ngeom,ndofn,pospg,wpg,N,Nxi,Neta,nenP,NP,model,perme_e) 
% [Ke,Ge,fe] = EleMat(Xe,ngeom,ndofn,pospg,wpg,N,Nxi,Neta,nenP,NP) 
% This function computes element matrices obtained for the 2D Stokes problem
% 
% Input:    
%   Xe:           Element nodes' coordinates
%   ngeom:        Number of nodes describing element geometry 
%                 (it is the number of element nodes except in PI+, which has 4 nodes and ngeom=3)
%   ndofn:        Number of degrees of freedom for velocity 
%   pospg, wpg:   Gauss points and weights int he reference element
%   N,Nxi,Neta:   Shape functions for velocity and its derivatives on Gauss points (local coordinates)
%   nenP:         Number of element nodes for the pressure field
%   NP:           Shape functions for the pressure field
%
% Output:
%   Ke: element viscosity matrix obtained by discretizing (grad v, grad u)_e
%   Ge: element gradient matrix (incompressibility constraint) obtained by discretizing (p,grad u)_e
%   fe: element sourde vector obtained by discretizing (v,f)_e
%
sourceTerm   =@model.source;
nu           = model.viscosity;
%permeability =@model.perme;

% Allocation
Ke = zeros(ndofn,ndofn); 
Ge = zeros(nenP,ndofn); 
fe = zeros(ndofn,1); 

% Number of Gauss points
ngaus = size(pospg,1); 

% Loop on Gauss points
for igaus = 1:ngaus
    % Shape functions on Gauss point igaus
    N_igaus    = N(igaus,:);  
    Nxi_igaus  = Nxi(igaus,:);   
    Neta_igaus = Neta(igaus,:);
    NP_ig      = NP(igaus,:);
    % Gauss point in global coordinates
    Xg = Isopar(Xe,N_igaus(1:ngeom));
    % Perme on gauss points
    %perme = permeability(Xg); %analytical fun
    perme = Isopar(perme_e,N_igaus(1:ngeom)); % evaluate nodal perme on gauss point
    % Jacobian matrix on the Gauss point
    Jacob = [Nxi_igaus(1:ngeom)*(Xe(1:ngeom,1))	    Nxi_igaus(1:ngeom)*(Xe(1:ngeom,2))   
             Neta_igaus(1:ngeom)*(Xe(1:ngeom,1))	Neta_igaus(1:ngeom)*(Xe(1:ngeom,2))];
    dvolu=wpg(igaus)*det(Jacob); 
    % Shape functions' derivatives in global coordinates
    res = Jacob\[Nxi_igaus;Neta_igaus]; 
    % Gradient
    Nx = [reshape([1;0]*res(1,:),1,ndofn); reshape([0;1]*res(1,:),1,ndofn)];
    Ny = [reshape([1;0]*res(2,:),1,ndofn); reshape([0;1]*res(2,:),1,ndofn)];
    % Divergence
    dN = reshape(res,1,ndofn);
    % Contribution to element matrices
    Ke = Ke + nu*(Nx'*Nx+Ny'*Ny)*dvolu; 
    Ng = [reshape([1;0]*N_igaus,1,ndofn); reshape([0;1]*N_igaus,1,ndofn)];
    if(perme>-eps)
        Me = nu/perme*(Ng'*Ng)*dvolu;
        Ke = Ke + Me;
    end

    Ge = Ge - NP_ig'*dN*dvolu; 
    
    f_igaus = sourceTerm(Xg); 
    Ngp = [reshape([1;0]*N_igaus,1,ndofn); reshape([0;1]*N_igaus,1,ndofn)];
    fe = fe + Ngp'*f_igaus*dvolu; 
end