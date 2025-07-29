function [mesh_groundAndBuildings]=generateSurface_buildings(mesh2D,fields,...
    minElevFromGround,translationVector,...
    doTreatSpecialBuildings,typeTreatSpecialBuildings,...
    meshFacadesWithQuads,idealizationLevel,...
    exportStepsToParaview)
fprintf('Generating surface mesh at building idealization level...\n')

verticalStepSize = minElevFromGround; % vertical extrusion steps to construct buildings

tic; fprintf('   Compute height field projections...\n')
[Z_ground]=projectMeshToGround(mesh2D,fields);
[heightRegion,heighRoofElement]=projectElementWiseFieldAppoximation(mesh2D,fields);
% warning('-----> he generat les alsades randomly <----------------------------------------')
% Z_ground = zeros(size(mesh2D.X,1),1);
% numRegions = mesh2D.numFields;
% heightRegion = rand(numRegions,1)*100;

elapsedTime = toc; fprintf('      ...%4.1f sec\n',elapsedTime);

X2D_ground = mesh2D.X;
Xground = [X2D_ground Z_ground];
groundRegion = mesh2D.groundRegion;


%% nou codi
buildingToBlock = mesh2D.buildingToBlock;

mesh_groundAndBuildings = mesh2D;
mesh_groundAndBuildings.X = Xground;
mesh_groundAndBuildings.T = mesh2D.T;

mesh_groundAndBuildings.facadeElementsRegion = cell(buildingToBlock.numBlocks,1);

% ===> per cada block <===
% global kkk;
% kkk = 0;
for iblock = (groundRegion+1):buildingToBlock.numBlocks
    
    mesh_groundAndBuildings = insertBuildingsInMesh(mesh_groundAndBuildings,...
        iblock,buildingToBlock,verticalStepSize,heightRegion,minElevFromGround);
    
end

if(isfield(mesh_groundAndBuildings,'facadeNodeData'))
    listNodes = unique(mesh_groundAndBuildings.T);
    oldNumNodes = size(mesh_groundAndBuildings.X,1);
    mesh_groundAndBuildings.X = mesh_groundAndBuildings.X(listNodes,:);
    newNumNodes = size(mesh_groundAndBuildings.X,1);
    mapOldToNew = zeros(oldNumNodes,1);
    mapOldToNew(listNodes) = 1:newNumNodes;
    mesh_groundAndBuildings.T = mapOldToNew(mesh_groundAndBuildings.T);


    mesh_groundAndBuildings.facadeNodeData.nodeList = ...
        mapOldToNew(mesh_groundAndBuildings.facadeNodeData.nodeList);
    mesh_groundAndBuildings.facadeNodeData.localNode_to_roofNode = ...
        mapOldToNew(mesh_groundAndBuildings.facadeNodeData.localNode_to_roofNode);
    mesh_groundAndBuildings.facadeNodeData.localNode_to_groundNode = ...
        mapOldToNew(mesh_groundAndBuildings.facadeNodeData.localNode_to_groundNode);

    % if the node does not have a ground node (it is between two roofs) then get the value of the roof
    listNotToGround = find(mesh_groundAndBuildings.facadeNodeData.localNode_to_groundNode==0);
    mesh_groundAndBuildings.facadeNodeData.localNode_to_groundNode(listNotToGround) = ...
        mesh_groundAndBuildings.facadeNodeData.localNode_to_roofNode(listNotToGround);
else
    fprintf('NO BUILDINGS/BLOCKS PRESENT IN MESH')
    error('no buildings')
end

%% things that are required frmo the code at some point (copied from generateSurface_block)

mesh_groundAndBuildings.terrain = zeros(length(mesh_groundAndBuildings.elementField),1);
groundElems = find(mesh_groundAndBuildings.elementField==mesh_groundAndBuildings.groundRegion);
mesh_groundAndBuildings.terrain(groundElems) = 1;
roofElements = find(mesh_groundAndBuildings.elementField>mesh_groundAndBuildings.groundRegion);
mesh_groundAndBuildings.terrain(roofElements) = 2;
facadeElems = find(mesh_groundAndBuildings.elementField==0);
mesh_groundAndBuildings.terrain(facadeElems) = 2;
maxZelems = max( [mesh_groundAndBuildings.X(mesh_groundAndBuildings.T(:,1),3) mesh_groundAndBuildings.X(mesh_groundAndBuildings.T(:,2),3) mesh_groundAndBuildings.X(mesh_groundAndBuildings.T(:,3),3)],[],2 );
tolSea = 0.1; %in meters
seaElems = find(maxZelems<tolSea);
mesh_groundAndBuildings.terrain(seaElems) = 0;


minAreaAllowedFacade = 0.0;
[validElements,areaElements] = checkMeshValidity(...
    mesh_groundAndBuildings,mesh_groundAndBuildings,...
    minElevFromGround,minAreaAllowedFacade,size(mesh2D.T,1));

facadeElements = size(mesh2D.T,1):size(mesh_groundAndBuildings.T,1);
invalidElements = facadeElements( find(~validElements( facadeElements ) ) );
numInvalidElements = length(invalidElements);

%mesh_groundAndBuildings.validity = validElements;
mesh_groundAndBuildings.area = areaElements;


if(exportStepsToParaview)
    mesh_groundAndBuildings.name = [mesh_groundAndBuildings.name '_facade3D'];
    Xsave = mesh_groundAndBuildings.X;
    mesh_groundAndBuildings.X = bsxfun(@minus,mesh_groundAndBuildings.X,[translationVector 0.0]);
    exportTriMeshToParaview(mesh_groundAndBuildings)
    mesh_groundAndBuildings.X = Xsave;
end


if(numInvalidElements>0)
    error('Invalid surface elements')
else    
    %mesh_groundAndBuildings = rmfield(mesh_groundAndBuildings,'validity');
    %fprintf('Closed surface mesh correctly defined...\n')
    
    mesh_groundAndBuildings.EN =giveConnectivity_ElementToNode(...
        mesh_groundAndBuildings.T,size(mesh_groundAndBuildings.X,1));
end


%%


% mesh_groundAndBuildings.name = 'provaTerraEdificis';
% exportTriMeshToParaview(mesh_groundAndBuildings)
% 
% if(exportStepsToParaview)
%     mesh_groundAndBuildings.name = [mesh_groundAndBuildings.name '_facade3D'];
%     Xsave = mesh_groundAndBuildings.X;
%     mesh_groundAndBuildings.X = bsxfun(@minus,mesh_groundAndBuildings.X,[translationVector 0.0]);
%     exportTriMeshToParaview(mesh_groundAndBuildings)
%     mesh_groundAndBuildings.X = Xsave;
% end


%% CUTRE TOT DE COP
% 
% 
% elemField = mesh2D.elementField;
% roofElements = find(elemField>groundRegion);
% numRoofElems = length(roofElements);
% heightRegion = heightRegion((groundRegion+1):end);
% steps = round(heightRegion./verticalStepSize);
% steps = unique(steps);
% steps = sort(steps);
% steps = setdiff(steps,0);
% minStep = min(steps);
% maxStep = max(steps);
% numLayers = maxStep-minStep+1;
% 
% heightRegionSteps = steps*verticalStepSize;
% 
% % disp('ei aixo esta falsejat')
% % heightRegionSteps = 10*(1:length(heightRegionSteps));
% % Z_ground = zeros(size(Z_ground));
% 
% numSrfElems = size(mesh2D.T,1);
% numVerticalElementsMax = numLayers-1;
% numPrismsMax = numVerticalElementsMax*numSrfElems;
% numNodesPrism = 6;
% TprismBig = zeros(numPrismsMax,numNodesPrism);
% numPrisms = 0;
% 
% %nodesRoofs = mesh.T(roofElements,:);
% %numRoofPoints = unique(nodesRoofs(:));
% numPointsSrf = size(X2D_ground,1); %mesh2D.X
% numPointsVolMax = numPointsSrf*numLayers;
% Zbig = zeros(numPointsVolMax,1);
% Zbig(1:numPointsSrf,:) = Z_ground;
% Xbig = zeros(numPointsVolMax,3);                        % XXX no cal aquesta matriu gorda (ara la faig per comoditat
% Xbig(1:numPointsSrf,:) = [X2D_ground Z_ground];         % XXX no cal
% 
% for istep = 1:length(steps)%maxStep%minStep:maxStep
%    
%     %heightStep = istep*verticalStepSize;
%     %heightStep = istep*verticalStepSize;
%     %heightStep = steps(istep)*verticalStepSize;
%     heightStep = heightRegionSteps(istep);
%     
%     nodes_previousLayer = (1:numPointsSrf)+ (numPointsSrf*(istep-1));
%     nodes_currentLayer = (1:numPointsSrf)+ (numPointsSrf*(istep));
%     
%     % extrudar tots els punts de superficie (potser no cal els que no
%     % necessitem, pero ocupar lespai per simplificar)
%     Zbig(nodes_currentLayer,:) = heightStep;
%     Xbig(nodes_currentLayer,1:2) = X2D_ground;          % XXX no cal
%     Xbig(nodes_currentLayer,3) = heightStep;            % XXX no cal
%     
%     % indexar els punts a partir dels de superficie per cada capa
%     % usar nomes elque toca per fer fer els prismes
%     for iregion = (groundRegion+1):mesh2D.numFields
%         if(heightRegion(iregion-1)>= heightStep - eps)
%             elementsInField = mesh2D.fieldElements{iregion};
%             for ielem=1:length(elementsInField)
%                 theElem = elementsInField(ielem);
%                 nodesElem = mesh2D.T(theElem,:);
% 
%                 numPrisms=numPrisms+1;
%                 TprismBig(numPrisms,1:3) = nodes_previousLayer(nodesElem);
%                 TprismBig(numPrisms,4:6) = nodes_currentLayer(nodesElem);
%             end
%         end
%     end
%     % eliminar despres els que no interessen de volume
% end


%%
% numBlocks           = mesh2D.buildingToBlock.numBlocks;
% mapBuildingToBlock  = mesh2D.buildingToBlock.map;
% blockElems          = mesh2D.buildingToBlock.blockElems;
% buildingsInBlock    = mesh2D.buildingToBlock.buildingsInBlock;
%%
% ZZ = setdiff(mesh_groundAndBuildings.facadeNodeData.nodeList,mesh_groundAndBuildings.T(:));
% if(not(isempty(ZZ)))
%     error('a')
% end
% ZZ = setdiff(mesh_groundAndBuildings.facadeNodeData.localNode_to_roofNode,mesh_groundAndBuildings.T(:));
% if(not(isempty(ZZ)))
%     error('b')
% end
% ZZ = setdiff(mesh_groundAndBuildings.facadeNodeData.localNode_to_groundNode,mesh_groundAndBuildings.T(:));
% if(not(isempty(ZZ)))
%     error('c')
% end


%%
% meshPri.T = TprismBig(1:numPrisms,:);
% meshPri.X = Xbig;
% meshPri.name = 'provaPrism';
% exportPrismMeshToParaview(meshPri)
% 
% meshSRF = mesh2D;
% meshSRF.X = [X2D_ground Z_ground];
% meshSRF.name = 'provaTerra';
% exportTriMeshToParaview(meshSRF)

