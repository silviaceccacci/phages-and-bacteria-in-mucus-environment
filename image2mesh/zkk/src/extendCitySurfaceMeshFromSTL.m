function [mesh_facade3D,domainLimits,domainLimitsIni]=extendCitySurfaceMeshFromSTL(...
    mesh_facade3D,domainLimits,edgeLength,elementOptions,options)

%% Extend
fprintf('Extending ground mesh...\n')
[mesh_facade3D,domainLimits]=extendMesh(mesh_facade3D,domainLimits,edgeLength,elementOptions,options);

domainLimitsIni = domainLimits;

%% Recompute ground marks
fprintf('Recomputing ground marks...\n')
[mesh_facade3D]=recomputeGroundMarks(mesh_facade3D);
warning('recomputeGroundMarks malament, tornar a posar')

%% Generate extra fields not generated for the STL case
mesh_facade3D.fieldElements = cell(2,1);
mesh_facade3D.fieldElements{mesh_facade3D.groundRegion}=find(mesh_facade3D.elementField==mesh_facade3D.groundRegion);
mesh_facade3D.fieldElements{mesh_facade3D.roofRegion}  =find(mesh_facade3D.elementField==mesh_facade3D.roofRegion);

mesh_facade3D.facadeElementsRegion = cell(2,1);
mesh_facade3D.facadeElementsRegion{mesh_facade3D.roofRegion} = find(mesh_facade3D.elementField==mesh_facade3D.facadeRegion);

mesh_facade3D.terrain = zeros(length(mesh_facade3D.elementField),1);
groundElems = find(mesh_facade3D.elementField==mesh_facade3D.groundRegion);
mesh_facade3D.terrain(groundElems) = 1;
roofElements = find(mesh_facade3D.elementField>mesh_facade3D.groundRegion);
mesh_facade3D.terrain(roofElements) = 2;
facadeElems = find(mesh_facade3D.elementField==0);
mesh_facade3D.terrain(facadeElems) = 2;
maxZelems = max( [mesh_facade3D.X(mesh_facade3D.T(:,1),3) mesh_facade3D.X(mesh_facade3D.T(:,2),3) mesh_facade3D.X(mesh_facade3D.T(:,3),3)],[],2 );
%tolSea = 0.1; %in meters
%seaElems = find(maxZelems<tolSea);
seaElems = [];
mesh_facade3D.terrain(seaElems) = 0;

mesh_facade3D.numFields = 2;

%%
mesh_facade3D.name = [mesh_facade3D.name '_facade3D'];
exportTriMeshToParaview(mesh_facade3D)

end


function [mesh,domainLimits]=extendMesh(mesh,domainLimits,edgeLength,elementOptions,options)

    disp('Set polybox from STL')
    tic
    [polyBox,domainLimits,zextend]=setPolyBoxSTL(mesh,domainLimits,edgeLength);
    toc
    
%     warning('fer que quan mallo el 2D no li passi tota la X sino nomes la dels boundaries i despres juntar, sino tira molts warnings..')
%     polyBox = rmfield(polyBox,'pointHole');
   

%     options.meshSize
%     options.meshAngle
    options = rmfield(options,'meshSize');
    options = rmfield(options,'meshAngle');
    options.meshAngle = 5;

    disp('Mesh surrounding area')
    options.edgeMarks = false;
    options.removeNotUsedNodes = false;
    planarExtendedMesh = meshConformingPolylines(polyBox,elementOptions,options);
    
%     numRepeatedNodes = 0;
%     if(isfield(polyBox,'mapToMesh'))
%         % not all nodes form the STL are in the box and we have to map them
%         numRepeatedNodes = length(polyBox.mapToMesh);
%         mapToMesh = [polyBox.mapToMesh
%                      1+(1:size(mesh.X,1))'];
%         planarExtendedMesh.T = mapToMesh(planarExtendedMesh.T);
%     end
    
%     planarExtendedMesh.name = [mesh.name '_planarExten'];
%     exportTriMeshToParaview(planarExtendedMesh)
    
%     figure;
%     plot(polyBox.X(polyBox.T(:),1),polyBox.X(polyBox.T(:),2),'*')
%     axis equal
%     hold on
%     plot(planarExtendedMesh.X(:,1),planarExtendedMesh.X(:,2),'ro')
%     plot(mesh.X(:,1),mesh.X(:,2),'go','MarkerSize',5)
%     title('boundary of stl and added boundary box')
        
    % Merge meshes
    numNodes_0 = size(mesh.X,1);
    numNodes_1 = size(planarExtendedMesh.X,1);
    Xnew = planarExtendedMesh.X;
%     mesh
%     size( zextend*ones(numNodes_1-numNodes_0,1))

    if(isfield(polyBox,'mapToMesh'))
        % not all nodes form the STL are in the box and we have to map them
        numRepeatedNodes = length(polyBox.mapToMesh);
        numNewNodesExtend = size(planarExtendedMesh.X,1)-numRepeatedNodes;
        mapToMesh = [polyBox.mapToMesh
                     size(mesh.X,1)+(1:numNewNodesExtend)'];
        planarExtendedMesh.T(:,:) = mapToMesh(planarExtendedMesh.T);
        
        Xnew = [mesh.X(:,1:2)  ; planarExtendedMesh.X( (numRepeatedNodes+1):end,:) ];
        Xnew = [Xnew  [ mesh.X(:,3) ; zextend*ones(numNewNodesExtend,1) ] ];
    else
        Xnew = [Xnew   [ mesh.X(:,3) ; zextend*ones(numNodes_1-numNodes_0,1) ] ];
    end
%     if(numRepeatedNodes==0)
%         Xnew = [Xnew   [ mesh.X(:,3) ; zextend*ones(numNodes_1-numNodes_0,1) ] ];
%     else
%         Xnew = [mesh.X  ; planarExtendedMesh.X( (numRepeatedNodes+1):end,:) ];
%         Xnew = [Xnew  [ mesh.X(:,3) ; zextend*ones(numNodes_1-numNodes_0-numRepeatedNodes,1) ] ];
%     end
    mesh.X = Xnew;
    
    
    mesh.extendedElementList = size(mesh.T,1) + (1:size(planarExtendedMesh.T,1));
    
    mesh.T = [mesh.T ; planarExtendedMesh.T];
    mesh.elementField = [mesh.elementField ; mesh.groundRegion*ones(size(planarExtendedMesh.T,1),1)];
    mesh.element = planarExtendedMesh.element;
    
    mesh.name = [mesh.name '_exten'];
    exportTriMeshToParaview(mesh)
    
%     mesh.fieldElements = [mesh.fieldElements ; mesh.groundRegion*ones(size(T,1),1)];
end

function [polyBox,domainLimits,zextend]=setPolyBoxSTL(mesh,domainLimits,edgeLength)

    doRemoveInterior = true;

    element.type = 'tri';
    element.numVertices = 3;
    element.numEdges = 3;
    element.order = 1;

    [boundNodes,boundElements,boundEdges] = giveBoundaryFromConnectivity(...
        mesh.T,size(mesh.X,1),element);

    zextend = min(mesh.X(boundNodes,3));
    
	numBoundElems = length(boundElements);
    TboundEdges = zeros(numBoundElems,2);
    for ibelem =1:numBoundElems
        theElem = boundElements(ibelem);
        theEdge = boundEdges(ibelem);
        TboundEdges(ibelem,:) = getEdge_tri(mesh.T(theElem,:),1,theEdge);
    end
    
    if(doRemoveInterior)
        numBoundNodes = length(boundNodes);
        mapToMesh = zeros(numBoundNodes,1);
        mapToMesh(:) = boundNodes;
        polyBox.mapToMesh = mapToMesh;
        
        mapMeshToPoly = zeros(size(mesh.X,1),1);
        mapMeshToPoly(boundNodes) = 1:numBoundNodes;
        TboundEdges(:,:) = mapMeshToPoly(TboundEdges);

        Xmesh = mesh.X(boundNodes,1:2);

        %warning('here doing ok the polybox... ;-)')
    else
        Xmesh = mesh.X(:,1:2);
    end
    
    spacingX = edgeLength;
    spacingY = edgeLength;

    xmax = max(mesh.X(:,1));
    xmin = min(mesh.X(:,1));
    ymax = max(mesh.X(:,2));
    ymin = min(mesh.X(:,2));
    
    if(isfield(domainLimits,'marginRight'))
        xmax = xmax+domainLimits.marginRight;
        ymax = ymax+domainLimits.marginUp;
        xmin = xmin-domainLimits.marginLeft;
        ymin = ymin-domainLimits.marginDown;
    else
        xmax = xmax+domainLimits.margin;
        ymax = ymax+domainLimits.margin;
        xmin = xmin-domainLimits.margin;
        ymin = ymin-domainLimits.margin;
    end
    
    nx = ceil((xmax-xmin)/spacingX);
    spacingX = (xmax-xmin)/nx - 1e-12;
    ny = ceil((ymax-ymin)/spacingY);
    spacingY = (ymax-ymin)/ny - 1e-12;
    
    Xpoints = xmin:spacingX:xmax; Xpoints = Xpoints';
    Ypoints = ymin:spacingY:ymax; Ypoints = Ypoints';
    
    mx = length(Xpoints)-1;
    my = length(Ypoints)-1;
    Xbox = [Xpoints(1:mx)              Ypoints(1)*ones(mx,1)
            Xpoints(mx+1)*ones(my,1)   Ypoints(1:my)
            Xpoints((mx+1):-1:2)       Ypoints(my+1)*ones(mx,1)
            Xpoints(1)*ones(my,1)      Ypoints((my+1):-1:2)];
        
% 	polyBox2D.X = Xbox;
         
	T = [ [1:size(Xbox,1)]' [2:size(Xbox,1) 1]'];
    
%     polyBox2D.T = T;
        
    T = T + size(Xmesh,1);
    
%     polyBox2D.globIdsXbox = (size(X,1)+1):(size(X,1)+size(Xbox,1));
    
    polyBox.X = [ Xmesh; Xbox];
    polyBox.T = [TboundEdges; T];
    polyBox.mark = zeros(size(polyBox.T,1),1);
    
%     buildingMark = 2;
%     boundaryMark = 3;
%     edgeMarks = [buildingMark*ones(size(polylines.T,1),1) ; boundaryMark*ones(size(T,1),1)];
%     
%     polyBox.buildingMark = buildingMark;
%     polyBox.boundaryMark = boundaryMark;
%     polyBox.mark = edgeMarks;
%     polyBox2D.buildingMark = buildingMark;
%     polyBox2D.boundaryMark = boundaryMark;
%     polyBox2D.mark = edgeMarks;
    
    pointHole = sum(mesh.X(:,1:2),1)/size(mesh.X,1);
    polyBox.pointHole = pointHole;
    
    domainLimits.type = 'reset';
    domainLimits.xmin = xmin;
    domainLimits.xmax = xmax;
    domainLimits.ymin = ymin;
    domainLimits.ymax = ymax;
    domainLimits.hbound = (spacingX+spacingY)/2.0;
end






%      roofElements = 1:size(mesh.T,1);
%      roofElements = roofElements(find(regions(roofElements)==mesh.roofRegion));
%      find(MAE(roofElements,1)
     

