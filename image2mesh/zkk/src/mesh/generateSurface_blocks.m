function [mesh_facade3D]=generateSurface_blocks(mesh2D,fields,...
    minElevFromGround,translationVector,...
    doTreatSpecialBuildings,typeTreatSpecialBuildings,...
    meshFacadesWithQuads,idealizationLevel,...
    exportStepsToParaview)
fprintf('Generating surface mesh at block idealization level...\n')

count = 1; maxIt = 10;

minAreaAllowedFacade = 0.0;

mesh2D.specialBuildings.doTreat  = doTreatSpecialBuildings;
mesh2D.specialBuildings.typeTreat= typeTreatSpecialBuildings;

lastNonFacadeElement = size(mesh2D.T,1);
%nonFacadeElements = 1:lastNonFacadeElement;
while( (count ==1 || numInvalidElements>0) && count<=maxIt )
    
    fprintf('   ...it %d\n',count)
    
    if(count>1)
        nonRoofRegions = unique(mesh_facade3D.facadeAssociatedRegion(invalidElements));
        nonRoofRegions = nonRoofRegions(find(nonRoofRegions>0));
        [mesh2D] = removeRegionsAsRoofs(mesh2D,nonRoofRegions) ;
        fprintf('      ...Num invalid removed regions: %d\n',length(nonRoofRegions))
    end
     
    tic; fprintf('   Generating facades...\n')
    mesh_facade = facadeDuplication(mesh2D,idealizationLevel);
    elapsedTime = toc; fprintf('      ...%4.1f sec\n',elapsedTime);

    facadeElements = (lastNonFacadeElement+1):size(mesh_facade.T,1);
 
    tic; fprintf('   Compute height field projections...\n')
    [mesh_facade3D]=projectingMeshToField('elementalStructured',mesh_facade,fields);
    elapsedTime = toc; fprintf('      ...%4.1f sec\n',elapsedTime);

    [validElements,areaElements] = checkMeshValidity(...
        mesh_facade3D,mesh_facade,minElevFromGround,minAreaAllowedFacade,size(mesh2D.T,1));
    
    mesh_facade3D.validity = validElements;
    mesh_facade3D.area = areaElements;
    
    invalidElements = facadeElements( find(~validElements( facadeElements ) ) );
    numInvalidElements = length(invalidElements);
    
    fprintf('      ...Num invalid elements: %d\n',numInvalidElements)
    count = count+1;
    
% mesh_facade3D.name =  [mesh_facade3D.name '_' int2str( count)];
% exportTriMeshToParaview(mesh_facade3D);
end

mesh_facade3D.terrain = zeros(length(mesh_facade3D.elementField),1);
groundElems = find(mesh_facade3D.elementField==mesh_facade3D.groundRegion);
mesh_facade3D.terrain(groundElems) = 1;
roofElements = find(mesh_facade3D.elementField>mesh_facade3D.groundRegion);
mesh_facade3D.terrain(roofElements) = 2;
facadeElems = find(mesh_facade3D.elementField==0);
mesh_facade3D.terrain(facadeElems) = 2;
maxZelems = max( [mesh_facade3D.X(mesh_facade3D.T(:,1),3) mesh_facade3D.X(mesh_facade3D.T(:,2),3) mesh_facade3D.X(mesh_facade3D.T(:,3),3)],[],2 );
tolSea = 0.1; %in meters
seaElems = find(maxZelems<tolSea);
mesh_facade3D.terrain(seaElems) = 0;


if(exportStepsToParaview)
    %exportTriMeshToParaview(mesh_facade3D);
    Xsave = mesh_facade3D.X;
    mesh_facade3D.X = bsxfun(@minus,mesh_facade3D.X,[translationVector 0.0]);
    exportTriMeshToParaview(mesh_facade3D)
    mesh_facade3D.X = Xsave;
end

[isConsistent,mesh_facade3D]=checkMeshConsistency_usedNodes(mesh_facade3D);
% if(checkMeshConsistency_usedNodes(mesh_facade3D))
%     error('Some lost nodes (not all nodes are present in the connectivity')
% end

if(numInvalidElements>0)
    error('Invalid surface elements')
else    
    mesh_facade3D = rmfield(mesh_facade3D,'validity');
    fprintf('Closed surface mesh correctly defined...\n')
    
    mesh_facade3D.EN =giveConnectivity_ElementToNode(mesh_facade3D.T,size(mesh_facade3D.X,1));
end

if(meshFacadesWithQuads)
    numFacTris = size(mesh_facade3D.T,1)-size(mesh2D.T,1);
    numFacQuads = numFacTris/2;
    numHybridElements = size(mesh2D.T,1)+numFacQuads;
    Thybrid = zeros(numHybridElements,4);
    elemNumNodes = zeros(numHybridElements,1);
    elemNumNodes(1:size(mesh2D.T,1)) = 3;
    Thybrid(1:size(mesh2D.T,1),1:3) = mesh_facade3D.T(1:size(mesh2D.T,1),:);
    elemNumNodes((size(mesh2D.T,1)+1):end) = 4;
    Thybrid( (size(mesh2D.T,1)+1):end ,1:4) = ...
        [mesh_facade3D.T((size(mesh2D.T,1)+1):2:end,1:2) mesh_facade3D.T((size(mesh2D.T,1)+2):2:end,2:3) ];
    mesh_facade3D.Thybrid = Thybrid;
    mesh_facade3D.hybrid = true;
end