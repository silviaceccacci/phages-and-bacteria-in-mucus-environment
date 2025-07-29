function [surfaceMesh]=buildClosedBoxSurfaceMesh(...
    topHeightOverTerrain,mesh_facade3D,elementOptions,parameters_tet,options)

topZ = max(mesh_facade3D.X(:,3));
ceilBoxZ = topZ+topHeightOverTerrain;

options.meshSize = parameters_tet.areaCeil;%10*parameters_tet.sizeCeil;
options.meshAngle = 15;
options.edgeMarks = false;

myEps = 1e-10;
xp = mesh_facade3D.X(:,1);
yp = mesh_facade3D.X(:,2);
minX = min(xp);
minY = min(yp);
maxX = max(xp);
maxY = max(yp);
n1 = find(abs(yp-minY)<myEps);% n1 = n1(find(yp(n1)<minY+myEps));
n2 = find(abs(xp-maxX)<myEps);% n1 = n1(find(yp(n1)<minY+myEps));
n3 = find(abs(yp-maxY)<myEps);
n4 = find(abs(xp-minX)<myEps);
[kk,i] = sort(xp(n1));  n1 = n1(i);
[kk,i] = sort(yp(n2));  n2 = n2(i);
[kk,i] = sort(-xp(n3)); n3 = n3(i);
[kk,i] = sort(-yp(n4)); n4 = n4(i);
n1 = n1(1:(end-1));
n2 = n2(1:(end-1));
n3 = n3(1:(end-1));
n4 = n4(1:(end-1));
polyBoundary.X = mesh_facade3D.X([n1; n2; n3; n4],1:2);
nn = size(polyBoundary.X,1);
polyBoundary.T = [[1:nn]' [2:nn 1]'];
polyBoundary.mark = zeros(size(polyBoundary.T,1),1);
[meshCeiling]=meshConformingPolylines(polyBoundary,elementOptions,options);
XCeil3D = [ meshCeiling.X ceilBoxZ*ones(size(meshCeiling.X,1),1)];

%downNodes = [polyBox2D.globIdsXbox'; polyBox2D.globIdsXbox(1)];
%upNodes   = size(mesh_facade3D.X,1) + ([1:size(polyBox2D.X,1) 1])';
downNodes = [n1; n2; n3; n4; n1(1)];
upNodes   = size(mesh_facade3D.X,1) + ([1:size(polyBoundary.X,1) 1])';

%Tlateral1 = [downNodes(2:end)  downNodes(1:(end-1))   upNodes(1:(end-1))  ];
Tlateral1 = [downNodes(2:end)  upNodes(1:(end-1))  downNodes(1:(end-1))   ];
Tlateral2 = [upNodes(2:end)    upNodes(1:(end-1))     downNodes(2:end)];
Tlateral = zeros(size(Tlateral1,1)+size(Tlateral2,1),3);
Tlateral(1:2:size(Tlateral,1),:) = Tlateral1;
Tlateral(2:2:size(Tlateral,1),:) = Tlateral2;

ceilElemId = -2;
lateralBoxFaces = -3;

surfaceMesh = mesh_facade3D;
surfaceMesh.X = [mesh_facade3D.X; XCeil3D];
surfaceMesh.T = [mesh_facade3D.T; meshCeiling.T+size(mesh_facade3D.X,1); Tlateral];%Tlateral1 ; Tlateral2];

%% Geometry facets
numFacetsRoofs = length(surfaceMesh.fieldElements) - surfaceMesh.groundRegion;
numFacetsFacades = length(surfaceMesh.facadeElementsRegion) - surfaceMesh.groundRegion;
numFacetsBound = 5;
numFacets = 1 + numFacetsRoofs + numFacetsFacades + numFacetsBound;

if(surfaceMesh.groundRegion~=1)
    error('not a problem but why is it changed?')
end
surfaceMesh.fieldElements{surfaceMesh.groundRegion} = ...
    find(surfaceMesh.elementField==surfaceMesh.groundRegion);
surfaceMesh.facets = cell(numFacets,1);
surfaceMesh.facetBoundId = zeros(numFacets,1);

surfaceMesh.facets(surfaceMesh.groundRegion:(numFacetsRoofs+1)) = ...
    surfaceMesh.fieldElements;%(surfaceMesh.groundRegion:end);
surfaceMesh.facetBoundId(surfaceMesh.groundRegion:(numFacetsRoofs+1)) = 5;

numFacetsIni = 1+numFacetsRoofs+1;
numFacetsEnd = numFacetsIni+numFacetsFacades-1;
surfaceMesh.facets(numFacetsIni:numFacetsEnd) = ...
    surfaceMesh.facadeElementsRegion((surfaceMesh.groundRegion+1):end);
surfaceMesh.facetBoundId(numFacetsIni:numFacetsEnd) = 5;
%% Alya Boundary markers for BC condition in 3D mesh
boundId = zeros(size(surfaceMesh.T,1),1);
% floor
numFloorElems = size(mesh_facade3D.T,1);
boundId(1:numFloorElems) = 5;
% ceiling
firstId = numFloorElems+1;
lastId = numFloorElems+size(meshCeiling.T,1);
boundId(firstId:lastId) = 6;
numFacetsIni = numFacetsEnd+1; numFacetsEnd = numFacetsIni; 
surfaceMesh.facets{numFacetsIni:numFacetsEnd} = (firstId:lastId)';
surfaceMesh.facetBoundId(numFacetsIni:numFacetsEnd) = 6;
%laterals
firstId = lastId+1;
lastId = lastId+length(n1)*2;
boundId(firstId:lastId) = 1;
numFacetsIni = numFacetsEnd+1; numFacetsEnd = numFacetsIni; 
surfaceMesh.facets{numFacetsIni:numFacetsEnd} = (firstId:lastId)';
surfaceMesh.facetBoundId(numFacetsIni:numFacetsEnd) = 1;
firstId = lastId+1;
lastId = lastId+length(n2)*2;
boundId(firstId:lastId) = 2;
numFacetsIni = numFacetsEnd+1; numFacetsEnd = numFacetsIni; 
surfaceMesh.facets{numFacetsIni:numFacetsEnd} = (firstId:lastId)';
surfaceMesh.facetBoundId(numFacetsIni:numFacetsEnd) = 2;
firstId = lastId+1;
lastId = lastId+length(n3)*2;
boundId(firstId:lastId) = 3;
numFacetsIni = numFacetsEnd+1; numFacetsEnd = numFacetsIni; 
surfaceMesh.facets{numFacetsIni:numFacetsEnd} = (firstId:lastId)';
surfaceMesh.facetBoundId(numFacetsIni:numFacetsEnd) = 3;
firstId = lastId+1;
lastId = lastId+length(n4)*2;
boundId(firstId:lastId) = 4;
numFacetsIni = numFacetsEnd+1; numFacetsEnd = numFacetsIni; 
surfaceMesh.facets{numFacetsIni:numFacetsEnd} = (firstId:lastId)';
surfaceMesh.facetBoundId(numFacetsIni:numFacetsEnd) = 4;

surfaceMesh.elementBoundId = boundId;


%% Sizing 
if(isfield(parameters_tet,'size') && parameters_tet.size)
    
    hnodes = zeros(size(surfaceMesh.X,1),1);
    
    if(isfield(parameters_tet,'sizeGround'))
        hnodes(1:size(mesh_facade3D.X,1)) = parameters_tet.sizeGround;
    else
        %disp('COMPUTE MIN EDGE LENGTH ADJACENT TO NODE NOT BELONGING TO FACADE')
        
        % Compensate the fact that if two sizes have been set, the mesh for
        % the coarse area is not generated taking into account the
        % factor_srf_volume. That is, the mesh is generated with the final
        % size and not mupliplied by factor_srf_volume (because otherwise
        % it is too coarse and the topography is not properply captured)
        areaVector = mesh_facade3D.area;
        areaVector(mesh_facade3D.coarseElements) = areaVector(mesh_facade3D.coarseElements)*...
            parameters_tet.factor_srf_volume^2;
        
        groundRegion = mesh_facade3D.groundRegion;
        for inode=1:size(mesh_facade3D.X,1)
            adjElems = find(mesh_facade3D.EN(:,inode));
            adjElemsRegion = mesh_facade3D.elementField(adjElems);
            adjElems_notGroundOrFac = adjElems(find(adjElemsRegion>=groundRegion));
            meanArea = sum(areaVector(adjElems_notGroundOrFac))/length(adjElems_notGroundOrFac);
            if(length(adjElems_notGroundOrFac)>0)
                
                edgeLength = sqrt(meanArea*2);
                hnodes(inode) = edgeLength/parameters_tet.factor_srf_volume;%(edgeLength^3)/(6.0*sqrt(2));
                
            else
                
                isFacadePoint = min(adjElemsRegion==mesh_facade3D.facadeRegion);
                
                if(isFacadePoint) 
                    % do nothing, we let tetgen choose the size
                    %hnodes(inode) = 5.0;
                    hnodes(inode) = -1.0;
                else
                    inode
                    length(unique(mesh_facade3D.T))
                    size(mesh_facade3D.X,1)
                    find(mesh_facade3D.EN(:,inode))
                    if(size(mesh_facade3D.X,1)~=length(unique(mesh_facade3D.T)))
                        error('aquiiiii')
                    end
                    find(mesh_facade3D.T(:)==inode)

                    meshSRF.X = mesh_facade3D.X;
                    meshSRF.T = mesh_facade3D.T(find(mesh_facade3D.EN(:,inode)),:);
                    figure(132); hold on;
                    for i=1:size(meshSRF.T,1)
                        xtri = meshSRF.X(meshSRF.T(i,:),:);
                        plot3(xtri(:,1),xtri(:,2),xtri(3,:),'*') 
                    end
                    listNodes = unique(meshSRF.T(:))
                    mapNodesKK = zeros(size(meshSRF.X,1),1);
                    mapNodesKK(listNodes) = 1:length(listNodes);
                    meshSRF.X = meshSRF.X(listNodes,:);
                    meshSRF.T = mapNodesKK(meshSRF.T);
                    meshSRF.name = ['error_mesh_facade3D'];
                    exportTriMeshToParaview(meshSRF)

                    meshSRF.X = mesh_facade3D.X;
                    meshSRF.T = mesh_facade3D.T;
                    meshSRF.name = ['error_mesh_facade3D_complete'];
                    exportTriMeshToParaview(meshSRF)

                    error('no adjacent elements')
                end
            end
        end
%         % for the floor nodes:
%         for inode=1:size(mesh_facade3D.X,1)
%             adjElems = find(mesh_facade3D.EN(:,inode));
%             adjElems = adjElems(find(mesh_facade3D.elementField(adjElems)>=groundRegion));
%             meanArea = sum(mesh_facade3D.area(adjElems))/length(adjElems);
%             edgeLength = sqrt(meanArea*2);
%             hnodes(inode) = edgeLength/parameters_tet.factor_srf_volume;%(edgeLength^3)/(6.0*sqrt(2));
%         end
%         % for the ceil nodes, imposed initially
        
        % interpolate size field on interior nodes on the facades
        if(isfield(mesh_facade3D,'facadeNodeData'))
            IF_nodes = mesh_facade3D.facadeNodeData.nodeList;
            num_IF = length(IF_nodes);
            IF_roof   = mesh_facade3D.facadeNodeData.localNode_to_roofNode;
            IF_ground = mesh_facade3D.facadeNodeData.localNode_to_groundNode;

            distToGround = mesh_facade3D.X(IF_nodes,3)-mesh_facade3D.X(IF_ground,3);
            distToRoof = mesh_facade3D.X(IF_roof,3)-mesh_facade3D.X(IF_nodes,3);

            distTotal = (distToGround + distToRoof);
            alpha_ground = distToGround./distTotal;
            beta_roof = distToRoof./distTotal;

            h_ground = hnodes(IF_ground);
            h_roof = hnodes(IF_roof);

            hvalue = alpha_ground.*h_ground + beta_roof.*h_roof;

            list_innerBlockInnerFacadeNodes = find(distTotal<eps);
            hvalue(list_innerBlockInnerFacadeNodes) = h_roof(list_innerBlockInnerFacadeNodes);

            if(find(isnan(hvalue)))
                alpha_ground
                h_ground 
                beta_roof
                h_roof
                error('not properly assigned size on facade node')
            end

            hnodes(IF_nodes) = hvalue;
        end
    end
    %hnodes(1:size(mesh_facade3D.X,1))
    
    hnodes((size(mesh_facade3D.X,1)+1):end) = parameters_tet.edgeLengthCeil;
    surfaceMesh.hnodes = hnodes;
    
    notAssignedNodes = find(hnodes<0);
    if(length(notAssignedNodes)>0)
        error('some size not assigned');
    end
end

%% Renumbering for output
surfaceMesh.elementField =...
   [surfaceMesh.elementField-surfaceMesh.groundRegion
    ceilElemId*ones(size(meshCeiling.T,1),1)
    lateralBoxFaces*ones(2*size(Tlateral1,1),1)];
surfaceMesh.groundRegion =0;
surfaceMesh.name = [mesh_facade3D.name '_withBox'];

sky = 3;
surfaceMesh.terrain =...
   [surfaceMesh.terrain
    sky*ones(size(meshCeiling.T,1),1)
    sky*ones(2*size(Tlateral1,1),1)];

% We will have:
% -  1:numFields for the roofs
% -  0 for the ground
% - -1 for the facades
% - -2 for the ceil
% - -3 for the lateral faces
surfaceMesh.numFields = surfaceMesh.numFields-1;

end





