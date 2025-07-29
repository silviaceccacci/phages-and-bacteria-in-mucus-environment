function Cv = CreConv_par_faster(velo,X,IEN,elem,isparall)
% Cv = CreConv(velo,X,IEN,elem);
% This function computes convection matrix
% for a two-dimensional Navier-Stokes problem
%
% Input:
%
%   velo:   velocity field
%   X:      nodal coordinates (velocity)
%   IEN:    connectivities (velocity)
%   elem:   type of element (0: quadrilateral, 1: triangle)
%
[N,Nxi,Neta,pospg,pespg,numel,nen,numnp,ngaus,nnodes,ndofn,nunk]=...
    getGaussPointsShapeFun_conv(elem,IEN,X);

% Allocates storage for matrices
Te=(1:ndofn)';
listi=zeros(ndofn*ndofn,1);%zeros(0,1);
listj=zeros(ndofn*ndofn,1);%zeros(0,1);
for iaux=1:ndofn
    %listi=[listi;Te(:)];
    %listj=[listj;Te(iaux)*ones(ndofn,1)];
    count_ini= 1+(ndofn*(iaux-1));
    count_end= 0+(ndofn*(iaux  ));
    listi(count_ini:count_end)=Te;
    listj(count_ini:count_end)=Te(iaux)*ones(ndofn,1);
end

CC = zeros(ndofn,ndofn,numel);
% Loop on elements
if(isparall)
    parfor ielem = 1:numel   
        % Xe: element nodes' coordinates 
        Xe = X(IEN(ielem,1:nnodes),:); 
        % Ve: element nodes' velocity
        Ve = velo(IEN(ielem,:),:); 
        % element matrix
        CC(:,:,ielem)=eleMat_conv(Xe,ndofn,N,Nxi,Neta,Ve,pespg,nnodes);
    end 
else
    for ielem = 1:numel   
        % Xe: element nodes' coordinates 
        Xe = X(IEN(ielem,1:nnodes),:); 
        % Ve: element nodes' velocity
        Ve = velo(IEN(ielem,:),:); 
        % element matrix
        CC(:,:,ielem)=eleMat_conv(Xe,ndofn,N,Nxi,Neta,Ve,pespg,nnodes);
    end 
end 

iiv = zeros(ndofn*ndofn*numel,1);
jjv = zeros(ndofn*ndofn*numel,1);
kkv = zeros(ndofn*ndofn*numel,1);
for ielem = 1:numel 
    % Te: current velocity element
    Te = reshape([IEN(ielem,:); numnp+IEN(ielem,:)],1,ndofn); 

    count_ini= 1+(ndofn*ndofn*(ielem-1));
    count_end= 0+(ndofn*ndofn*(ielem  ));
    iiv(count_ini:count_end) = Te(listi);
    jjv(count_ini:count_end) = Te(listj);
    kkv(count_ini:count_end) = CC(:,:,ielem);
end
clear CC

% Assembly
Cv = sparse(iiv,jjv,kkv,nunk,nunk);    
    
end
