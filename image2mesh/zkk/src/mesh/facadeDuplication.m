function [mesh] = facadeDuplication(mesh,idealizationLevel)

    if(strcmp(idealizationLevel,'building'))
        error('At building level, the srf mesh generation is by constructing the volume buildings')
        [mesh] = facadeDuplication_fromBoundaryEdges(mesh);
    else
        [mesh] = facadeDuplication_fromBoundaryRoofs(mesh);
        %[mesh] = facadeDuplication_fromBoundaryEdges_illes(mesh);
        %[mesh] = facadeDuplication_fromConnectivity(mesh);
    end
end

function [mesh] = facadeDuplication_fromBoundaryRoofs(mesh)
%% get boundary of roof regions
elementRegion = mesh.elementField;
roofElements = find(elementRegion>mesh.groundRegion);

facadeElementsRegion = cell(length(mesh.fieldElements),1);

[matrixAdjacentElement,matrixLocalFaceAdjacentElement,boundaryNodes]=...
    getMatrixAdjacentElement_general_city(mesh.T(roofElements,:),size(mesh.X,1),mesh.element);

boundElem = find( min(matrixAdjacentElement,[],2) == 0);
numBoundElem = length(boundElem);

maxNumBuildingEdges = 2*numBoundElem;
Tfacades1 = zeros(maxNumBuildingEdges,3); % elements adjacen to ground
Tfacades2 = zeros(maxNumBuildingEdges,3); % elements adjacen to roof

numNodes = size(mesh.X,1);
mapOldToNew = (1:numNodes)';
isNodeToDuplicate = false(numNodes,1);

Xdupli = zeros(maxNumBuildingEdges,2);

facadeAssociatedRegion = zeros(maxNumBuildingEdges,1);

numIniElements = size(mesh.T,1);
countFacadeElements = 0;
countNewNodes = numNodes;
locduplid = 0;
elemEdges = [ 1 2; 2 3; 3 1];
for ielem=1:length(boundElem)
    iroofelem = boundElem(ielem);
    globelem = roofElements(iroofelem);
    for iedge = 1:size(matrixAdjacentElement,2)
        if(matrixAdjacentElement(iroofelem,iedge)==0)
            
            bnodes = mesh.T(globelem,elemEdges(iedge,:));
            if(isNodeToDuplicate(bnodes(1))==false)
                isNodeToDuplicate(bnodes(1)) = true;
                countNewNodes = countNewNodes+1;
                mapOldToNew(bnodes(1)) = countNewNodes;
                
                locduplid = locduplid+1;
                Xdupli(locduplid,:) = mesh.X(bnodes(1),:);
            end
            if(isNodeToDuplicate(bnodes(2))==false)
                isNodeToDuplicate(bnodes(2)) = true;
                countNewNodes = countNewNodes+1;
                mapOldToNew(bnodes(2)) = countNewNodes;
                
                locduplid = locduplid+1;
                Xdupli(locduplid,:) = mesh.X(bnodes(2),:);
            end
            
            countFacadeElements = countFacadeElements+1;
            Tfacades1(countFacadeElements,:) = ...
                 [bnodes([1]) mapOldToNew(bnodes(1)) bnodes([2])];
                 %[bnodes([2 1]) mapOldToNew(bnodes(1))];
                 %[bnodes([1 2]) mapOldToNew(bnodes(1))];
            Tfacades2(countFacadeElements,:) = ...
                 [bnodes(2) mapOldToNew(bnodes(2)) mapOldToNew(bnodes(1))];
             
            facadeAssociatedRegion(countFacadeElements) = elementRegion(globelem);
            
            facadeElementsRegion{elementRegion(globelem)} =...
                [facadeElementsRegion{elementRegion(globelem)};...
                (numIniElements+countFacadeElements)
                (numIniElements+countFacadeElements+1)];
        end
    end
end

Tfacades = zeros(2*countFacadeElements,3);
Tfacades(1:2:size(Tfacades,1),:) = Tfacades1(1:countFacadeElements,:);
Tfacades(2:2:size(Tfacades,1),:) = Tfacades2(1:countFacadeElements,:);

facadeAssociatedRegionAll = zeros(2*countFacadeElements,1);
facadeAssociatedRegionAll(1:2:size(Tfacades,1)) = facadeAssociatedRegion(1:countFacadeElements);
facadeAssociatedRegionAll(2:2:size(Tfacades,1)) = facadeAssociatedRegion(1:countFacadeElements);

Tnew= mesh.T;
Tnew(roofElements,:) = mapOldToNew(Tnew(roofElements,:));
newT = [Tnew ; Tfacades];

totNewNodes =locduplid;% countNewNodes-numNodes;
Xdupli = Xdupli(1:totNewNodes,:);
newX = [ mesh.X ; Xdupli];

facadeAssociatedRegion = [zeros(size(Tnew,1),1); facadeAssociatedRegionAll];

%%

facadeRegion = mesh.groundRegion - 1;
newElementRegions = facadeRegion*ones(size(Tfacades,1),1);
mesh.facadeRegion = facadeRegion;

mesh.facadeElementsRegion = facadeElementsRegion;

mesh.T = newT;
mesh.X = newX;
mesh.elementField = [mesh.elementField ; newElementRegions];

mesh.X = [mesh.X zeros(size(mesh.X,1),1)];

mesh.facadeAssociatedRegion = facadeAssociatedRegion;

mesh.X(mesh.T(roofElements,:),3) = 50;

mesh.name = [mesh.name '_facade'];

end

%% wrong or slow funcitons
function [mesh] = facadeDuplication_fromBoundaryEdges_illes(mesh)
%% get info from edges and mesh
edgeInfo = mesh.edges;
Tedges = edgeInfo.T;
NNedgeId = edgeInfo.NNedgeId;
edgeMark = edgeInfo.mark;

ENmesh = mesh.EN;
elementRegion = mesh.elementField;
roofElements = find(elementRegion>mesh.groundRegion);
%% get building edges and update them (some buildings may have been removed from catastral info)

buildingEdges = find(edgeMark== edgeInfo.buildingMark);
numIniEdges = length(buildingEdges);

disp('mirar de vectoritzar aquesta part!!!!')

fprintf('   ...%d initial forced building edges',numIniEdges)

edgeMark_updated = zeros(numIniEdges,1);
for iedge = 1:numIniEdges
    n1 = Tedges(iedge,1);
    n2 = Tedges(iedge,2);
    adjElems = find( ENmesh(:,n1) );
    adjElems = adjElems( find( ENmesh(adjElems,n2) ) );
    if(max(elementRegion(adjElems))==mesh.groundRegion)
        edgeMark_updated(iedge) = edgeInfo.buildingMark;
    end
end

buildingEdges = find(edgeMark_updated == edgeInfo.buildingMark);
numBuildingEdges = length(buildingEdges);
globalEdgeIdToBuilding = zeros(size(Tedges,1),1);
globalEdgeIdToBuilding(buildingEdges) = 1:numBuildingEdges;

fprintf(' reduced to %d.\n',numBuildingEdges)


%% duplicate boundary nodes
TedgesBuildings = Tedges(buildingEdges,:);
nodesToDuplicate = unique(TedgesBuildings(:));
numNodesToDuplicate = length(nodesToDuplicate);

numNodes = size(mesh.X,1);
mapOldToNew = (1:numNodes)';
mapOldToNew(nodesToDuplicate) = (numNodes+1):(numNodes+numNodesToDuplicate);

Tfacades1 = zeros(numBuildingEdges,3); % elements adjacen to ground
Tfacades1(:,[1 2 3]) = ...
    [ TedgesBuildings(:,[1 2]) mapOldToNew(TedgesBuildings(:,1))];
Tfacades2 = zeros(numBuildingEdges,3); % elements adjacen to roof
Tfacades2(:,[1 2 3]) = ...
    [ TedgesBuildings(:,2) mapOldToNew(TedgesBuildings(:,2)) mapOldToNew(TedgesBuildings(:,1))];

% % reorientate facades if inverted 
% elemEdges = [ 1 2; 2 3; 3 1];
% for ielem = 1:length(roofElements)
%     theElem = roofElements(ielem);
%     for iedge =1:3
%         edgeNodes = mesh.T(theElem,elemEdges(iedge,:));
%         edgeId = NNedgeId(edgeNodes(1),edgeNodes(2));
%         if(edgeId>0)
%             globalEdgeIdToBuilding(edgeId)
%             buildEdgeId = globalEdgeIdToBuilding(edgeId);
%             if(TedgesBuildings(buildEdgeId,1)==edgeNodes(1)) %if the edge is given in the roof orientation, then invert
%                 Tfacades1(buildEdgeId,:) = Tfacades1(buildEdgeId,[2 1 3]);
%                 Tfacades2(buildEdgeId,:) = Tfacades2(buildEdgeId,[1 3 2]);
%             end
%         end
%     end
% end

Tfacades = zeros(2*numBuildingEdges,3);
Tfacades(1:2:size(Tfacades,1),:) = Tfacades1;
Tfacades(2:2:size(Tfacades,1),:) = Tfacades2;

%error('AIXO NO TEN EN COMPTE Q HI HA DUPLICACIO A ROOOOOFFFFFSSSS!!!!')
Tini = mesh.T;
Tnew = Tini;
%Troof = Tini(roofElements,:);
%Troof(:) = mapOldToNew(Troof(:));
%Tini(roofElements,:) = Troof;
for iregion=(mesh.groundRegion+1):mesh.numFields
    elementRegion = mesh.fieldElements{iregion};
    Tnew(elementRegion,:) = reshape(mapOldToNew(Tini(elementRegion,:)),length(elementRegion),3);
    keepOldIndex = mapOldToNew(Tini(elementRegion,:));
    mapOldToNew(keepOldIndex(:)) = keepOldIndex(:);
end

Xini = mesh.X;
Xdupli = Xini(nodesToDuplicate,:);
newX = [Xini; Xdupli];

newT = [Tini ; Tfacades];

%% Ensure no inverted facade faces
ENnew =giveConnectivity_ElementToNode(newT,size(newX,1));

%check for each pair of facade elements if they are properly oriented with
%respect to their adjacent roof element, and if not, invert their order
elemEdges = [ 1 2; 2 3; 3 1];
for ielem=1:2:size(Tfacades,1)
    ielemg = ielem;
    ielemr = ielem+1;
    elemg_glob = size(Tini,1) + ielemg;
    elemr_glob = size(Tini,1) + ielemr;
    
    bnodesRoof = Tfacades(ielemr,[2 3]);
    %find adjacent eleme to ielemr
    adjElem = ENnew(:,bnodesRoof(1));
    ENnew(adjElem,bnodesRoof(2))
    adjElem = adjElem(ENnew(adjElem,bnodesRoof(2)));
    adjElem = setdiff(adjElem,elemr_glob);
    inverted = true;
    for iedge = 1:size(elemEdges,2)
        edgeNodes = newT(adjElem,elemEdges(iedge,:));
        if(edgeNodes(2) == bnodesRoof(1) && edgeNodes(2) == bnodesRoof(1))
            inverted = false;
            break;
        end
    end
    if(inverted)
        newT(elemg_glob,:) = newT(elemg_glob,[2 1 3]);
        newT(elemr_glob,:) = newT(elemr_glob,[1 3 2]);
    end
end
mesh.EN = ENnew;

%%

facadeRegion = mesh.groundRegion - 1;
newElementRegions = facadeRegion*ones(size(Tfacades,1),1);
mesh.facadeRegion = facadeRegion;

mesh.T = newT;
mesh.X = newX;
mesh.elementField = [mesh.elementField ; newElementRegions];

mesh.X = [mesh.X zeros(size(mesh.X,1),1)];


mesh.X(mesh.T(roofElements,:),3) = 50;

mesh.name = [mesh.name '_facade'];

end


function [mesh] = facadeDuplication_fromConnectivity(mesh)
%% find boundary between regions
groundRegion = mesh.groundRegion;

maxNumBoundaries = 10*mesh.numFields;
maxNumBoundaryElementsRegion = 1000;

neighboringRegions = zeros(mesh.numFields,mesh.numFields);
numBoundaries = 0;
elementsOnBoundaries = zeros(maxNumBoundaries,maxNumBoundaryElementsRegion,2,2);
countBoundary = zeros(maxNumBoundaries,1);
for ielem=1:size(mesh.T,1)
    currentRegion = mesh.elementField(ielem);
    if(currentRegion>groundRegion)
        for iedge = 1:mesh.element.numEdges
            neighElem = mesh.matrixAdjacentElement(ielem,iedge);
            neighEdge = mesh.matrixLocalFaceAdjacentElement(ielem,iedge);
            if(neighElem >0)
                neighRegion = mesh.elementField(neighElem);
                if(neighRegion~=currentRegion)
                    if(neighboringRegions(currentRegion,neighRegion)==0)
                       numBoundaries = numBoundaries+1; 
                       neighboringRegions(currentRegion,neighRegion) = numBoundaries;
                       neighboringRegions(neighRegion,currentRegion) = numBoundaries;
                    end
                    boundId = neighboringRegions(currentRegion,neighRegion);
                    countBoundary(boundId) = countBoundary(boundId) + 1;

                    if(currentRegion<neighRegion)
                        elementsOnBoundaries(boundId,countBoundary(boundId),1,:) = [ ielem neighElem];
                        elementsOnBoundaries(boundId,countBoundary(boundId),2,:) = [ iedge neighEdge];
                    else
                        elementsOnBoundaries(boundId,countBoundary(boundId),1,:) = [ neighElem ielem];
                        elementsOnBoundaries(boundId,countBoundary(boundId),2,:) = [ neighEdge iedge];
                    end
                end
            end
        end
    end
end
countBoundary = countBoundary(1:numBoundaries);
maxNumBoundaryElementsRegion = max(countBoundary);
elementsOnBoundaries = elementsOnBoundaries(1:numBoundaries,1:maxNumBoundaryElementsRegion,1:2,1:2);
%% duplicate boundary nodes
% the region with the lowest id will keep the initial nodes, the other will
% have the duplicated nodes
currentNode = size(mesh.X,1);
oldNumNodes = currentNode;
totalNewNodes = sum(countBoundary);
newX = zeros(totalNewNodes,2);
X = mesh.X;
newX = [X; newX];
totalNewElems = totalNewNodes*2;
newT = zeros(totalNewElems,3);
T = mesh.T;
%newT = [T;newT];
newElements = 0 ;
mapOldNodeToNew = 1:size(mesh.X,1);
groundBoundaries = zeros(numBoundaries*2,1);
newNodeBound = 0;
for iboun = 1:numBoundaries
    newNodeList = [ (currentNode+1):(currentNode+countBoundary(iboun))   ,   currentNode+1];
    for iBedge = 1:countBoundary(iboun)
        elem2 = elementsOnBoundaries(iboun,iBedge,1,2);
        edge2 = elementsOnBoundaries(iboun,iBedge,2,2);
        
        edgeNodes2 = getEdge(T(elem2,:),mesh.element,edge2);
        
        mapOldNodeToNew(edgeNodes2(1)) =  newNodeList(iBedge);
        newX( newNodeList(iBedge),:) = X(edgeNodes2(1),:);
        
        newNodeBound = newNodeBound+1;
        groundBoundaries(newNodeBound) = edgeNodes2(1);
    end
    for iBedge = 1:countBoundary(iboun)
        elem1 = elementsOnBoundaries(iboun,iBedge,1,1);
        elem2 = elementsOnBoundaries(iboun,iBedge,1,2);
        edge1 = elementsOnBoundaries(iboun,iBedge,2,1);
        edge2 = elementsOnBoundaries(iboun,iBedge,2,2);
        
        edgeNodes1 = getEdge(T(elem1,:),mesh.element,edge1);
        edgeNodes2 = getEdge(T(elem2,:),mesh.element,edge2);
        
 
        newElements = newElements+1;
        newT(newElements,:) = [edgeNodes1([2 1]) mapOldNodeToNew(edgeNodes2(1))];
        newElements = newElements+1;
        newT(newElements,:) = [mapOldNodeToNew(edgeNodes2([2 1])) edgeNodes1(1)];
    end
    currentNode = currentNode+countBoundary(iboun);
end
newElementRegions = zeros(newElements,1);%-ones(newElements,1);

targetElems = find(mesh.elementField>groundRegion);
Ttarget = T(targetElems,:);
Ttarget(:) = mapOldNodeToNew(Ttarget(:));
T(targetElems,:) = Ttarget;

boundaryNodes.roof = oldNumNodes:(oldNumNodes+totalNewNodes);
boundaryNodes.ground = unique(groundBoundaries);
boundaryNodes.exterior = mesh.boundaryNodes;

mesh.boundaryNodes = boundaryNodes;

%newT(:) = mapOldNodeToNew(newT(:));

% for iboun = 1:numBoundaries
%     for iBedge = 1:countBoundary(iboun)
%         elem1 = elementsOnBoundaries(iboun,iBedge,1,1);
%         elem2 = elementsOnBoundaries(iboun,iBedge,1,2);
%         edge1 = elementsOnBoundaries(iboun,iBedge,2,1);
%         edge2 = elementsOnBoundaries(iboun,iBedge,2,2);
%         
%         edgeNodes1 = getEdge(T(elem1,:),mesh.element,edge1);
%         edgeNodes2 = getEdge(T(elem2,:),mesh.element,edge2);
%         
%         newElements = newElements+1;
%         newT(newElements,:) = [edgeNodes1([2 1]) edgeNodes2(1)];
%         
%         newElements = newElements+1;
%         newT(newElements,:) = [edgeNodes1([2 1]) edgeNodes2(1)];
%         %newT(newElements,:) = [edgeNodes2([2 1]) edgeNodes1(1)];
%     end
% end


newT = [T;newT];
%newT = T;

mesh.T = newT;%[mesh.T ; newT];
mesh.X = newX;%[mesh.X ; newX];
mesh.elementField = [mesh.elementField ; newElementRegions];

mesh.X = [mesh.X zeros(size(mesh.X,1),1)];
% for ielem = 1:size(mesh.T,1)
%     if(mesh.elementField>1)
%         mesh.X(mesh.T(ielem,:),3) = 50;
%     end
% end

mesh.X(mesh.T(targetElems,:),3) = 50;

mesh.name = [mesh.name '_facade'];

end
   