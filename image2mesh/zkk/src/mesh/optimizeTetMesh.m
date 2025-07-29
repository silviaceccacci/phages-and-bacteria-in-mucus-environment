function [mesh] = optimizeTetMesh(mesh,parameters)
    curveMeshPath = '/Users/abel/Dropbox/PHD/codi/codesPHD_matlab/codes/';
    %%
    addpath([curveMeshPath 'qualityMeshes_firstOrder_tets/'])
    addpath([curveMeshPath 'minimizationMethods/'])
    addpath([curveMeshPath 'numericalDerivationGeneral/'])
    addpath([curveMeshPath 'minimizationMethods/'])
    addpath([curveMeshPath 'qualityMeshes_firstOrder_linear/parametricObjectiveFunction/'])
    addpath([curveMeshPath 'qualityMeshes_firstOrder_linear/colourAdjacencyClassifyingFunctions/'])
%     addpath([curveMeshPath 'auxiliarSmoothing/'])
%     addpath([curveMeshPath 'auxiliarMeshes/'])
%     addpath([curveMeshPath 'qualityMeshes_highOrder/colourAdjacencyClassifying_HO/'])
    global withImages; withImages = false;
    %%
    minQ = 0.2%0.05;%0.1;
    nlevels = 2;
    maxIt = 10;
    %%
    if(~isfield(mesh,'boundaryNodes'))
        mesh.boundaryNodes = unique(mesh.boundaries.faces(:));
    end
    
    if(~isfield(mesh,'EN'))
        EN =giveConnectivity_ElementToNode(mesh.T,size(mesh.X,1));
    else
        EN = mesh.EN;
    end
    
    mesh.X = mesh.X';
    mesh.T = mesh.T(:,[1 3 2 4]);
    %%
    Winv = [ 1.0    -(1/3)*sqrt(3)       -(1/3)*sqrt(6)/2
             0.0     (2/3)*sqrt(3)       -(1/3)*sqrt(6)/2 
             0.0      0.0                        sqrt(6)/2  ];
	mesh.idealElements = Winv;
    %% Summary of optimization for linear tets
    [optimization,computeQuality]=setSmoothingProperties_tets();
    
    % % Compute quality
    fprintf('   Computing quality');
    quality = feval(computeQuality,mesh.X,mesh.T,Winv);
    fprintf(': min(%1.2f), mean(%1.2f), std(%1.2f),max(%1.2f)\n',min(quality),mean(quality),std(quality),max(quality));
    
    % % Get low quality elements and nodes to smooth    
    badElements = find(quality<minQ);
    badNodes = unique(mesh.T(badElements,:));
    badNodes = setdiff(badNodes,mesh.boundaryNodes);
    % % Add extra levels of nodes to give freedom to the optimization
    for ilevel=1:nlevels
        indexes = find(EN(:,badNodes));
        [connectedElems kk] = ind2sub(size(EN),indexes);
        badNodes = unique(mesh.T(unique(connectedElems),:));
        badNodes = setdiff(badNodes,mesh.boundaryNodes);
    end
    
    %mesh.nodesSmooth  = setdiff(badNodes,mesh.boundaryNodes);
    mesh.nodesSmooth  = badNodes;
    fprintf('   Num nodes to smooth: %d\n',length(mesh.nodesSmooth));
    
    % % Optimize mesh
    warning('Set a new smoothing framework, and not the preconditioner from HO')
    if(~isempty(mesh.nodesSmooth))
        [mesh.X,it,time] = feval(optimization,mesh,[],maxIt,Winv); % mainSmoothingGS_tets  
    end
    %% Calling older matlab code (slow)
%     options.subdomain = true;
%     %     
%     [mesh,mesh0] = smooth_3DLinear(mesh,options);
    %%
    fprintf('   Computing quality');
    mesh.quality = feval(computeQuality,mesh.X,mesh.T,Winv);
    fprintf(': min(%1.2f), mean(%1.2f), std(%1.2f),max(%1.2f)\n',min(quality),mean(quality),std(quality),max(quality));
    mesh.X = mesh.X';
    mesh.T = mesh.T(:,[1 3 2 4]);
end  
    
%%
    %[Winv_subElements] = giveWinvSubElements(element, linearMesh );
    %linearMeshToSmooth.T = linearMesh.T(:,[1,3,2,4]);