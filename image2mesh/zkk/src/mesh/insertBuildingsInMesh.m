function [mesh] = insertBuildingsInMesh(...
        mesh,iblock,buildingToBlock,verticalStepSize,heightRegion,minElevFromGround)
    
    minDistanceBetweenRoofs = minElevFromGround;
    
    % ===> agafo els buildings <===
    buildings = buildingToBlock.buildingsInBlock{iblock};
    elemsInBlock = [];
    for ibuild = 1:length(buildings)
        theBuilding = buildings(ibuild);
        if(theBuilding>size(mesh.fieldElements,1))
           theBuilding
            size(mesh.fieldElements)
        end
        elemsInBlock =[elemsInBlock; mesh.fieldElements{theBuilding}];
    end
    if(length(unique(elemsInBlock))<length(elemsInBlock))
        error('element with two marks')
    end   
    if( ne( length(elemsInBlock),length(buildingToBlock.blockElems{iblock}) ) )
        length(elemsInBlock)
        length(buildingToBlock.blockElems{iblock})
        error('The building to block structure has differed from the current marks, should be updated')
    end
        
    %elemsInBlock = buildingToBlock.blockElems{iblock};
    nodesBlockSrf = unique(mesh.T(elemsInBlock,:));
    Xsrf_block = mesh.X(nodesBlockSrf,:);
    maxZGround = max(Xsrf_block(:,3));
  
    %fprintf('\n ---> maxZGroudn %f\n',maxZGround)
    %fprintf('\n ---> length(elemsInBlock) %d\n',length(elemsInBlock))
    
    %% Collect buildings and heights (validate and order them)
    
    % ===> miro les seves alçades <===
    heightBuildings = heightRegion(buildings);

    % ===> les reestructuro en capes <===   
    [heightBuildings_sorted,sortedIndex] = sort(heightBuildings);
    
    numHeights = length(heightBuildings_sorted);
    isValidHeight = true(numHeights,1);
    %if( max(heightBuildings_sorted)<= maxZGround+minElevFromGround + eps )
    if( heightBuildings_sorted(end)<= maxZGround+minElevFromGround + eps )
        return % EXIT, NO BUILDING/BLOCK TO BE INSERTED
    else % the maximum height is already accepted
        for iheight = (numHeights-1):-1:1
            
            height1 = heightBuildings_sorted(iheight);
            height2 = heightBuildings_sorted(iheight+1);
                        
%             localBuildingToRemove = sortedIndex(iheight);
%             globalBuildingToRemove = buildings(localBuildingToRemove);
%             elemsInBuilding = mesh.fieldElements{globalBuildingToRemove};
%             nodesInBuilding = unique(mesh.T(elemsInBuilding,:));
%             maxZGround = max(mesh.X(nodesInBuilding,3));
            
            if(height1 <= maxZGround+minElevFromGround + eps) % check if above ground
                %disp('below ground')
                isValidHeight(iheight) = false;
                heightBuildings_sorted(iheight) = -1;%
                localBuildingToRemove = sortedIndex(iheight);
                globalBuildingToRemove = buildings(localBuildingToRemove);
                heightBuildings(localBuildingToRemove) = -1;
                %buildings = setdiff(buildings,buildingToRemove);
                elemsInBuilding = mesh.fieldElements{globalBuildingToRemove};
                elemsInBlock = setdiff(elemsInBlock,elemsInBuilding);
                nodesBlockSrf = unique(mesh.T(elemsInBlock,:));
                Xsrf_block = mesh.X(nodesBlockSrf,:);
                maxZGround = max(Xsrf_block(:,3)); 
                
                mesh.elementField(elemsInBuilding) = mesh.groundRegion; %<===============================removed elements from region
                
            elseif( height2-height1< minDistanceBetweenRoofs ) %check if similar to previous height
                %disp('above ground but similar to previous')
                isValidHeight(iheight) = false;
                
                heightBuildings_sorted(iheight) = height2;%(height1+height2)/2.0;
                heightBuildings_sorted(iheight+1) = height2;%(height1+height2)/2.0;
                
                localBuildingToRemove = sortedIndex(iheight);
                heightBuildings(localBuildingToRemove) = height2;%-1;
            end

        end
    end

    %fprintf('\n ---> length(elemsInBlock) %d\n',length(elemsInBlock))
    
    heightRegionSteps = heightBuildings_sorted(find(isValidHeight));
    numLayers = length(heightRegionSteps) +1;
    
    %% ===> les extrudo i construeixo els prismes <===

    numNodesBlockSrf = length(nodesBlockSrf); 
    
    numNodesMaxVolBlock = numNodesBlockSrf*numLayers;
    Xbig = zeros(numNodesMaxVolBlock,3);                        
    Xbig(1:numNodesBlockSrf,:) = Xsrf_block;

    numElemsBlock = length(elemsInBlock);
    numVerticalElementsMax = numLayers-1;
    numPrismsMax = numVerticalElementsMax*numElemsBlock;
    numNodesPrism = 6;
    TBuildingBig = zeros(numPrismsMax,numNodesPrism);
    numPrisms = 0;

    %%mapGlobalNodeToLocal(nodesBlockSrf) = 1:numNodesBlockSrf;
    mapGlobalNodeToLocal = sparse(nodesBlockSrf,ones(size(nodesBlockSrf)),...
        1:numNodesBlockSrf,max(nodesBlockSrf),1);

    numBuildings = length(buildings);
        
    %facadeNoes_list         = zeros(numNodesMaxVolBlock,1);
    %facadeNodes_toGround    = 1:numNodesBlockSrf;
    facadeNodes_toRoof      = 1:numNodesBlockSrf;
        
    for istep = 1:length(heightRegionSteps)

        heightStep = heightRegionSteps(istep);

        nodes_previousLayer = (1:numNodesBlockSrf) + (numNodesBlockSrf*(istep-1));
        nodes_currentLayer  = (1:numNodesBlockSrf) + (numNodesBlockSrf*(istep  ));

        % extrudar tots els punts de superficie (potser no cal els que no
        % necessitem, pero ocupar lespai per simplificar)
        %Zbig(nodes_currentLayer,:) = heightStep;
        Xbig(nodes_currentLayer,1:2) = Xsrf_block(:,1:2);%X2D_ground;          % XXX no cal
        Xbig(nodes_currentLayer,3) = heightStep;            % XXX no cal

        % indexar els punts a partir dels de superficie per cada capa
        % usar nomes elque toca per fer fer els prismes
        for ibuild_local = 1:numBuildings%(groundRegion+1):mesh.numFields    
            
            iregion = buildings(ibuild_local);

            if(heightBuildings(ibuild_local) >= heightStep - eps)

                elementsInField = mesh.fieldElements{iregion};

                for ielem=1:length(elementsInField)
                    theElem = elementsInField(ielem);
                    nodesElem = mesh.T(theElem,:);
                    nodesElem = mapGlobalNodeToLocal(nodesElem);

                    numPrisms=numPrisms+1;
                    TBuildingBig(numPrisms,1:3) = nodes_previousLayer(nodesElem);
                    TBuildingBig(numPrisms,4:6) = nodes_currentLayer(nodesElem);
                    
                    % anar actualitzant toRoof
                    facadeNodes_toRoof(nodesElem) = nodes_currentLayer(nodesElem);
                    
                end
            end
        end

        %usedNodes = unique(TprismBig);
        %mapLocalNodeToGlobal(usedNodes) = 

%         global kkk;
%         kkk = kkk+1
%         global XKK;
%         global TKK;
%         Taux = TBuildingBig(1:numPrisms,:);
%         TKK2 = [TKK; Taux + size(XKK,1)];
%         XKK2 = [XKK; Xbig];
%         meshPri.T = TKK2;
%         meshPri.X = XKK2;
% 
%     %     meshPri.T = TBuildingBig;
%     %     meshPri.X = Xbig;
%         textNum = int2str(kkk);
%         textNumLarge = '00000';
%         textNumLarge((end-length(textNum)+1):end) = textNum;
%         meshPri.name = ['prova/provaPrism' textNumLarge ];
%         exportPrismMeshToParaview(meshPri)
    end

    TBuildingBig = TBuildingBig(1:numPrisms,:);
            
%     heightBuildings
%    
%     global kkk;
%     kkk = kkk+1;
%     meshPri.T = TBuildingBig;
%     meshPri.X = Xbig;
% 
%     textNum = int2str(kkk);
%     textNumLarge = '00000';
%     textNumLarge((end-length(textNum)+1):end) = textNum;
%     meshPri.name = ['prova/provaPrism' textNumLarge ]
%     exportPrismMeshToParaview(meshPri)
    
    
    % eliminar despres els que no interessen de volume i de srf inicial
    element.type = 'pri';
    [matrixAdjacentElement,matrixLocalFaceAdjacentElement] =...
        getMatrixAdjacentElement_general(TBuildingBig,size(Xbig,1),element);
    [boundElements,boundFaces] = find(matrixAdjacentElement == 0);
    faceVertices = [1 3 2 0 
                4 5 6 0
                1 2 5 4  
                2 3 6 5
                3 1 4 6 ];
    numVerticesFaces = [3 3 4 4 4];
    
    numBElems = length(boundElements);
    %Tsurf = zeros(numBElems*2,3);
    Troof = zeros(numElemsBlock,3);
    Tfacade = zeros(numBElems*2,3);
    %countBound = 0;
    countRoof = 0;
    countFacade = 0;
    for ibound = 1:numBElems
        theElem = boundElements(ibound);
        iFace = boundFaces(ibound);
        if(iFace>1) % we dont accept the faces on the ground
            numVertices = numVerticesFaces(iFace);
            localNodes = faceVertices(iFace,1:numVertices);
            faceNodes = TBuildingBig(theElem,localNodes);

            if(numVertices==3)
                countRoof = countRoof + 1;
                Troof(countRoof,:) = faceNodes(1:3);
            elseif(numVertices==4)
                countFacade = countFacade + 1;
                Tfacade(countFacade,:) = faceNodes(1:3);
                countFacade = countFacade + 1;
                Tfacade(countFacade,:) = faceNodes([1 3 4]);
            end
            
%             countBound = countBound + 1;
%             Tsurf(countBound,:) = faceNodes(1:3);
%             if(numVertices==4)
%                 countBound = countBound + 1;
%                 Tsurf(countBound,:) = faceNodes([1 3 4]);
%             end
        end
    end
    Troof = Troof(1:countRoof,:);
    Tfacade = Tfacade(1:countFacade,:);
%     Tsurf = Tsurf(1:countBound,:);
      
    %validNodes = unique(Tsurf);
    validNodes = unique([Tfacade; Troof]);
    %numValidNodes = length(validNodes);
    %mapToValid = zeros(size(Xbig,1),1);
    %mapToValid(validNodes) = 1:numValidNodes;
    %Xvalid = 
    nodesToAdd = validNodes(find(validNodes>numNodesBlockSrf));
    numNodesToAdd = length(nodesToAdd);
    
    mapLocalNodeToGlobal = zeros(numNodesMaxVolBlock,1);
    mapLocalNodeToGlobal(1:numNodesBlockSrf) = nodesBlockSrf;
    mapLocalNodeToGlobal(nodesToAdd) = size(mesh.X,1) + (1:numNodesToAdd);
    
    mesh.X = [mesh.X
              Xbig(nodesToAdd,:)];
    %TsurfGlobal = mapLocalNodeToGlobal(Tsurf);
    Troof_global = mapLocalNodeToGlobal(Troof);
    Tfacade_global = mapLocalNodeToGlobal(Tfacade);
    
    mesh.T(elemsInBlock,:) = Troof_global;
    numPrevElems = size(mesh.T,1);
    mesh.T = [mesh.T; Tfacade_global];
    
    %mesh.elementField(elemsInBlock) = elemsBuildings_saved(elemsInBlock);
    mesh.elementField(elemsInBlock) = iblock;  % <==================================== no he actualitzat el field_elements
    
    facadeRegion = mesh.groundRegion - 1;
    %newElementRegions = facadeRegion*ones(size(Tfacades,1),1);
    mesh.facadeRegion = facadeRegion;
    mesh.elementField = [mesh.elementField; facadeRegion*ones(countFacade,1)];
    
    mesh.facadeElementsRegion{iblock} = numPrevElems:size(mesh.T,1);
    
    %%
    if(isfield(mesh,'facadeNodeData'))
        facadeNodeData = mesh.facadeNodeData;
    else
        facadeNodeData.nodeList = [];
        facadeNodeData.localNode_to_roofNode = [];
        facadeNodeData.localNode_to_groundNode = [];
    end
    %interiorFacadeNodes = setdiff(nodesToAdd,unique(Troof(:)));  
    facadeNodes = unique(Tfacade);
    facadeNodes_notGround = facadeNodes(find(facadeNodes>numNodesBlockSrf));
    interiorFacadeNodes = setdiff(facadeNodes_notGround,unique(Troof));
    interiorFacadeNodes_global = mapLocalNodeToGlobal(interiorFacadeNodes);
    groundNodes = mod(interiorFacadeNodes,numNodesBlockSrf);
    groundNodes(find(groundNodes==0)) = numNodesBlockSrf;
    groundNodes_global = mapLocalNodeToGlobal(groundNodes);
    roofNodes = facadeNodes_toRoof(groundNodes);
    roofNodes_global = mapLocalNodeToGlobal(roofNodes);
    
    facadeNodeData.nodeList = ...
        [facadeNodeData.nodeList
         interiorFacadeNodes_global];
    facadeNodeData.localNode_to_roofNode = ...
        [facadeNodeData.localNode_to_roofNode 
         roofNodes_global];
    facadeNodeData.localNode_to_groundNode = ...
        [facadeNodeData.localNode_to_groundNode 
         groundNodes_global];
    
     mesh.facadeNodeData = facadeNodeData;
     
    %%
    
%     meshSRF.T = Tsurf;
%     meshSRF.X = Xbig;
%     
%     global kkk;
%     kkk = kkk+1;
%     textNum = int2str(kkk);
%     textNumLarge = '00000';
%     textNumLarge((end-length(textNum)+1):end) = textNum;
%     meshSRF.name = ['prova/provaSRF' textNumLarge];
%     exportTriMeshToParaview(meshSRF)
    
% %     global kkk;
% %     kkk = kkk+1;
%     global XKK;
%     global TKK;
%     TKK = [TKK; TBuildingBig + size(XKK,1)];
%     XKK = [XKK; Xbig];
%     meshPri.T = TKK;
%     meshPri.X = XKK;
%     
% %     meshPri.T = TBuildingBig;
% %     meshPri.X = Xbig;
% 
% %     textNum = int2str(kkk);
% %     textNumLarge = '00000';
% %     textNumLarge((end-length(textNum)+1):end) = textNum;
% %     meshPri.name = ['prova/provaPrism' ];
% % 
% %    meshPri.name = ['prova/provaPrism' textNumLarge];
% %    exportPrismMeshToParaview(meshPri)
%     
%     %error('mira')
    
    
    
    %%
            
%         global kkk;
%         kkk = kkk+1;
%         global XKK; global TKK;
%         Taux =  TBuildingBig(1:numPrisms,:)+size(XKK,1);
%         TKK = [TKK; Taux];
%         XKK = [XKK ; Xbig];
%         meshPriKK.T = TKK;
%         meshPriKK.X = XKK;
%         meshPriKK.name = ['prova/provaPrism' int2str(kkk)];
%         exportPrismMeshToParaview(meshPriKK)
        
    %%
    %         for iheight = 1:(numHeights-1)
% 
%             height1 = heightBuildings_sorted(iheight);
%             height2 = heightBuildings_sorted(iheight+1);
%             if(height1 <= maxZGround+minElevFromGround + eps) 
%                 isValidHeight(iheight) = false;
%                 heightBuildings_sorted(iheight) = -1;%
%                 localBuildingToRemove = sortedIndex(iheight);
%                 globalBuildingToRemove = buildings(localBuildingToRemove);
%                 heightBuildings(localBuildingToRemove) = -1;
%                 %buildings = setdiff(buildings,buildingToRemove);
%                 elemsInBuilding = mesh.fieldElements{globalBuildingToRemove};
%                 elemsInBlock = setdiff(elemsInBlock,elemsInBuilding);
%                 nodesBlockSrf = unique(mesh.T(elemsInBlock,:));
%                 Xsrf_block = mesh.X(nodesBlockSrf,:);
%                 maxZGround = max(Xsrf_block,3);   
%             elseif( height2-height1< minDistanceBetweenRoofs )
%                 isValidHeight(iheight) = false;
%                 
%                 heightBuildings_sorted(iheight) = -1;%(height1+height2)/2.0;
%                 heightBuildings_sorted(iheight+1) = (height1+height2)/2.0;
%                 
%                 localBuildingToRemove = sortedIndex(iheight);
%                 heightBuildings(localBuildingToRemove) = -1;
%             end
% 
%         end
    %% ===> ho acoplo amb la X global <===
    
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
%     for iregion = (groundRegion+1):mesh_groundAndBuildings.numFields
%         if(heightRegion(iregion-1)>= heightStep - eps)
%             elementsInField = mesh_groundAndBuildings.fieldElements{iregion};
%             for ielem=1:length(elementsInField)
%                 theElem = elementsInField(ielem);
%                 nodesElem = mesh_groundAndBuildings.T(theElem,:);
% 
%                 numPrisms=numPrisms+1;
%                 TprismBig(numPrisms,1:3) = nodes_previousLayer(nodesElem);
%                 TprismBig(numPrisms,4:6) = nodes_currentLayer(nodesElem);
%             end
%         end
%     end
%     % eliminar despres els que no interessen de volume i de srf inicial