function [K,G,f] = CreMat_par_faster(X,IEN,XP,IENP,elemType,model,isParall) 
% [K,G,f] = CreMat(X,IEN,XP,IENP,elem) 
% This function computes viscosity matrix K, gradient matrix G
% and vector f (source term) obtained by discretizing the Galerkin
% weak form of the Stokes problem in 2D using finite elements.
%
% Input:
%
%   X:      nodal coordinates (velocity)
%   IEN:    connectivities (velocity)
%   XP:     nodal coordinates (pressure)
%   IENP:   connectivities (pressure)
%   elem:   type of element (0: quadrilateral, 1: triangle)

[N,Nxi,Neta,NP,NxiP,NetaP,pospg,wpg,numel,nen,nenP,numnp,ngaus,ndofn,nunk,nunkP,ngeom]=...
    getGaussPointsShapeFun(elemType,IEN,X,IENP,XP);

% Allocates storage for matrices
iiv = zeros(ndofn*ndofn*numel,1);
jjv = zeros(ndofn*ndofn*numel,1);
kkv = zeros(ndofn*ndofn*numel,1);
iip = zeros(nenP*ndofn*numel,1);
jjp = zeros(nenP*ndofn*numel,1);
kkp = zeros(nenP*ndofn*numel,1);

f = zeros(nunk,1); 

Te=(1:ndofn)';
listi=zeros(ndofn*ndofn,1);
listj=zeros(ndofn*ndofn,1);
for iaux=1:ndofn
    count_ini= 1+(ndofn*(iaux-1));
    count_end= 0+(ndofn*(iaux  ));
    listi(count_ini:count_end)=Te(:);
    listj(count_ini:count_end)=Te(iaux)*ones(ndofn,1);
end
Tep=(1:nenP)';
listip=zeros(0,1);
listjp=zeros(0,1);
for iaux=1:ndofn
    count_ini= 1+(nenP*(iaux-1));
    count_end= 0+(nenP*(iaux  ));
    listip(count_ini:count_end)=Tep;
    listjp(count_ini:count_end)=Te(iaux)*ones(nenP,1);
end

KK = zeros(ndofn,ndofn,numel);
GG = zeros(nenP,ndofn,numel);
ff = zeros(ndofn,numel);

if(isParall)
    parfor ielem = 1:numel 
        Xe = X(IEN(ielem,1:ngeom),:); 
        perme_e = model.perme(IEN(ielem,1:ngeom));
        [KK(:,:,ielem),GG(:,:,ielem),ff(:,ielem)] =...
            EleMat(Xe,ngeom,ndofn,pospg,wpg,N,Nxi,Neta,nenP,NP,...
            model,perme_e); 
    end
else
    for ielem = 1:numel 
        Xe = X(IEN(ielem,1:ngeom),:); 
        perme_e = model.perme(IEN(ielem,1:ngeom));
        [KK(:,:,ielem),GG(:,:,ielem),ff(:,ielem)] =...
            EleMat(Xe,ngeom,ndofn,pospg,wpg,N,Nxi,Neta,nenP,NP,...
            model,perme_e); 
    end
end

map_to_master_v = model.periodic_data.map_to_master_v;
map_to_master_p = model.periodic_data.map_to_master_p;
for ielem = 1:numel 
    elem_nodes = IEN(ielem,:);
    elem_nodes = map_to_master_v(elem_nodes);
    Te = reshape([elem_nodes; numnp+elem_nodes],1,ndofn); 

    TeP = IENP(ielem,:);
    TeP = map_to_master_p(TeP);

    count_ini= 1+(ndofn*ndofn*(ielem-1));
    count_end= 0+(ndofn*ndofn*(ielem  ));
    iiv(count_ini:count_end) = Te(listi);
    jjv(count_ini:count_end) = Te(listj);
    kkv(count_ini:count_end) = KK(:,:,ielem);

    count_ini= 1+(nenP*ndofn*(ielem-1));
    count_end= 0+(nenP*ndofn*(ielem  ));
    iip(count_ini:count_end) = TeP(listip);
    jjp(count_ini:count_end) = Te(listjp);
    kkp(count_ini:count_end) = GG(:,:,ielem);

    f(Te) = f(Te) + ff(:,ielem); 
end
clear KK; clear GG; clear ff;

K = sparse(iiv,jjv,kkv,nunk,nunk);
G = sparse(iip,jjp,kkp,nunkP,nunk);
