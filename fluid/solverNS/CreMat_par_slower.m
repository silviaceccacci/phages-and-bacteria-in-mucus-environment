function [K,G,f] = CreMat_par_faster(X,IEN,XP,IENP,elemType,model,periodic_data) 
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

% Preallocate storage for sparse matrices and force vector
iiv = zeros(ndofn*ndofn*numel, 1);
jjv = zeros(ndofn*ndofn*numel, 1);
kkv = zeros(ndofn*ndofn*numel, 1);
iip = zeros(nenP*ndofn*numel, 1);
jjp = zeros(nenP*ndofn*numel, 1);
kkp = zeros(nenP*ndofn*numel, 1);
f_i = zeros(ndofn*numel, 1); % Row indices for f
f_v = zeros(ndofn*numel, 1); % Values for f

% Define element connectivity indices
Te = (1:ndofn)';
listi = zeros(ndofn*ndofn, 1);
listj = zeros(ndofn*ndofn, 1);

for iaux = 1:ndofn
    count_ini = 1 + (ndofn * (iaux - 1));
    count_end = ndofn * iaux;
    listi(count_ini:count_end) = Te(:);
    listj(count_ini:count_end) = Te(iaux) * ones(ndofn, 1);
end

Tep = (1:nenP)';
listip = zeros(nenP * ndofn, 1);
listjp = zeros(nenP * ndofn, 1);

for iaux = 1:ndofn
    count_ini = 1 + (nenP * (iaux - 1));
    count_end = nenP * iaux;
    listip(count_ini:count_end) = Tep;
    listjp(count_ini:count_end) = Te(iaux) * ones(nenP, 1);
end

% Initialize cell arrays for thread-safe parallel assembly
iiv_cell = cell(numel, 1);
jjv_cell = cell(numel, 1);
kkv_cell = cell(numel, 1);
iip_cell = cell(numel, 1);
jjp_cell = cell(numel, 1);
kkp_cell = cell(numel, 1);
f_i_cell = cell(numel, 1);
f_v_cell = cell(numel, 1);

% Element matrices and force vector computation
map_to_master = periodic_data.map_to_master;
parfor ielem = 1:numel
    % Element geometry and permeability
    Xe = X(IEN(ielem, 1:ngeom), :);
    perme_e = model.perme(IEN(ielem, 1:ngeom));
    
    % Compute element stiffness, coupling, and force matrices
    [KK_e, GG_e, ff_e] = EleMat(Xe, ngeom, ndofn, pospg, wpg, N, Nxi, Neta, nenP, NP, model, perme_e);
    
    % Get global node indices
    elem_nodes = IEN(ielem, :);
    elem_nodes = map_to_master(elem_nodes);
    Te = reshape([elem_nodes; numnp + elem_nodes], 1, ndofn);

    TeP = IENP(ielem, :);
    TeP = map_to_master(TeP);

    % Assemble sparse matrix indices and values
    iiv_cell{ielem} = Te(listi);
    jjv_cell{ielem} = Te(listj);
    kkv_cell{ielem} = KK_e(:);

    iip_cell{ielem} = TeP(listip);
    jjp_cell{ielem} = Te(listjp);
    kkp_cell{ielem} = GG_e(:);

    % Assemble force vector indices and values
    f_i_cell{ielem} = Te(:);
    f_v_cell{ielem} = ff_e(:);
end

% Concatenate results from all workers
iiv = vertcat(iiv_cell{:});
jjv = vertcat(jjv_cell{:});
kkv = vertcat(kkv_cell{:});
iip = vertcat(iip_cell{:});
jjp = vertcat(jjp_cell{:});
kkp = vertcat(kkp_cell{:});
f_i = vertcat(f_i_cell{:});
f_v = vertcat(f_v_cell{:});

% Assemble global sparse matrices and force vector
K = sparse(iiv, jjv, kkv, nunk, nunk);
G = sparse(iip, jjp, kkp, nunkP, nunk);
f = accumarray(f_i, f_v, [nunk, 1]);

