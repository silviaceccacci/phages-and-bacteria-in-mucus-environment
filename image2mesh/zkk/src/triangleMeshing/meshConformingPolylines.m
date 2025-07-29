function [mesh]=meshConformingPolylines(polylines,elementOptions,options)

    if(~isfield(options,'removeNotUsedNodes'))
        removeNotUsedNodes = true;
    else
        removeNotUsedNodes = options.removeNotUsedNodes;
    end

    if(nargin==3 && isfield(options,'edgeMarks'))      
        readEdgeMarks = options.edgeMarks;
    else
        readEdgeMarks = true;
    end
    if(nargin==3 && isfield(options,'meshSize'))
        meshSize = options.meshSize;
        sizeString =[' -a' int2str(meshSize)];
        if(meshSize==0)
            sizeString = [];
        end
    else
        meshSize = 1000;
        sizeString =[' -a' int2str(meshSize)];
        % or sizeString = '';
    end
    
    if(isfield(polylines,'sizeBox'))
        sizeString = ' -a ';%[]; %it will be assigned for each region
    end
    
    if(nargin==3 && isfield(options,'meshAngle'))
        meshAngle = options.meshAngle;
        angleString =[' -q' int2str(meshAngle)];
        if(meshAngle==0)
            angleString = [];
        end
    else
        angleString =' -q10 ' ;
    end
    if(readEdgeMarks)
       stringEdge = ' -e ' ; 
    else
       stringEdge = ' ' ; 
    end
    
    if(nargin==4 && isfield(options,'constrainToBoundary'))
        if(options.constrainToBoundary)
            boundaryConstrinSwitch = '-Y';
        else
            boundaryConstrinSwitch = ' ';
        end
    else
        boundaryConstrinSwitch = ' ';
    end
    
    coarseMarkElements = 1000000;
    
    tempName = [options.tempFolder 'ztri_temp'];    
    [inpName] = buildTriangleInput(polylines,tempName,meshSize,coarseMarkElements);   
    logName = [tempName '.log'];
    
    triangleExe = '!./src/triangleMeshing/triangle/triangle.x ';
 
    optionsTriangle = [' -p ' inpName ' -A ' angleString  sizeString boundaryConstrinSwitch stringEdge ];
    
%     triangleMeshing = [ triangleExe optionsTriangle  ];%' &> ' logName ] ;
    triangleMeshing = [ triangleExe optionsTriangle  ' &> ' logName ] ;
%     fprintf('Output of triangle delaunay mesher: %s\n',logName)

    eval(triangleMeshing)
    
    mesh = readTriangleOutput( tempName, polylines, readEdgeMarks, coarseMarkElements );
    
    element.dim = 2;
    element.numCoord = 2;
    element.type = 'tri';
    element.order = 1;
    element.distribution = 'lineal';
    mesh.element = defineElement(element);
    
    mesh.element.coord = [ -1 -1;  1 -1 ; -1 1];
    
    mesh.element.quadrature = elementOptions.quadrature;
    mesh.element.orderQuadrature = elementOptions.orderQuadrature;
    mesh.element.shapeFunctions = elementOptions.shapeFunctions;
    
    mesh.EN =giveConnectivity_ElementToNode(mesh.T,size(mesh.X,1));

    if(removeNotUsedNodes)
        if(~checkMeshConsistency_usedNodes(mesh))
            usedNodes = unique(mesh.T);
            mapToNew = zeros(size(mesh.X,1),1);
            mapToNew(usedNodes) = 1:length(usedNodes);

            mesh.X = mesh.X(usedNodes,:);
            mesh.T = mapToNew(mesh.T); 
        end
    end
end


function [inpName] = buildTriangleInput(polylines,tempName,meshSize,coarseMarkElements)
    X = polylines.X;
    T = polylines.T;

    numNodeToDel = size(X,1);
    numEdges = size(T,1); % la idea es q aixo acabi sent una estructura per guardar regions i tal
    if(~isfield(polylines,'pointHole'))    
        numHoles = 0;
    else
        numHoles = 1;
        pointHole = polylines.pointHole;
    end
    inpName = [tempName '.poly'];
    stream = fopen(inpName,'w');
    
    fprintf(stream,'%d 2 0 0 0\n',numNodeToDel);
    fprintf(stream,'%d %f %f\n', [1:numNodeToDel; X']);
    fprintf(stream,'%d 1\n',numEdges);
    fprintf(stream,'%d %d %d %d\n', [1:numEdges ; T'; polylines.mark']);
    fprintf(stream,'\n');
    fprintf(stream,'%d\n',numHoles);
    if(numHoles>0)
        fprintf(stream,'%d %f %f\n', [1; pointHole']);
    end
    fprintf(stream,'\n');
    
    if(~isfield(polylines,'regions'))
        if(~isfield(polylines,'sizeBox'))
            numRegions = 1;
            atributeGroundRegion = 1;
            sizeRegion = -1;
            fprintf(stream,'%d\n',numRegions);
            fprintf(stream,'%d %f %f %f %f %f\n',[1 min(X(:,1))+0.1 min(X(:,2))+0.1  atributeGroundRegion sizeRegion]  );
        else
            if(~isfield(polylines,'extraRegions'))
                atributeGroundRegion = 1;
                fprintf(stream,'%d\n',2);
                fprintf(stream,'%d %f %f %d %f \n',[1 polylines.coordIn  atributeGroundRegion meshSize]  );
                fprintf(stream,'%d %f %f %d %f \n',[2 polylines.coordOut coarseMarkElements   polylines.sizeOut]  );
            else
                extraRegions = polylines.extraRegions;
                numExtraRegions = length(extraRegions.mark);
                extraMark = extraRegions.mark;
                extraCM = extraRegions.position;
                numRegions = 2+numExtraRegions;
                atributeGroundRegion = 1;
                fprintf(stream,'%d\n',numRegions);
                fprintf(stream,'%d %f %f %d %f \n',[1 polylines.coordIn  atributeGroundRegion meshSize]  );
                fprintf(stream,'%d %f %f %d %f \n',[2 polylines.coordOut coarseMarkElements   polylines.sizeOut]  );
                fprintf(stream,'%d %f %f %d %f \n',[3:numRegions
                                                    extraCM(:,1)'
                                                    extraCM(:,2)'
                                                    extraMark'%((3:numRegions)-1)%extraMark' 
                                                    meshSize*ones(1,numExtraRegions)    ]  );
            end
        end
    else
        numSpecialRegions = length(polylines.regions);
        numRegions = 1+numSpecialRegions;
        
        atributeGroundRegion = 1;
        sizeRegion = -1;
        fprintf(stream,'%d\n',numRegions);
        fprintf(stream,'%d %f %f %f %f %f\n',[1 min(X(:,1))+0.1 min(X(:,2))+0.1  atributeGroundRegion sizeRegion]  );
        
        atributeGroundRegion = 1; %we'll consider it as ground but we will sample it with lidar data
        for i=1:numSpecialRegions
            xcm = polylines.regions{i}.xcm;
            ycm = polylines.regions{i}.ycm;
            if(isfield(polylines.regions{i},'msize'))
                sizeRegion = polylines.regions{i}.msize;
            else
                sizeRegion = meshSize/10;
            end
            fprintf(stream,'%d %f %f %f %f %f\n',[(1+i) xcm ycm  atributeGroundRegion sizeRegion]  );
        end
    end
    %Optional line: <# of regional attributes and/or area constraints>
    %Optional following lines: <region #> <x> <y> <attribute> <maximum area>
    
% 	polylines.numRegions = iregion-1;
%     polylines.segmentRegion = segmentRegion;
%     polylines.numRegionSegments = numRegionSegments;
%     polylines.regionSegments = regionSegments;
    
%     numRegions = polylines.numRegions;
%     atributeRegionStreet = 0;
%     sizeRegion = -1;
%     fprintf(stream,'%d\n',numRegions);
%     %fprintf(stream,'%d %f %f %f %f %f\n',[1 min(X(:,1))+0.1 min(X(:,2))+0.1  atributeRegionStreet sizeRegion]  );
%     for i=1:polylines.numRegions
%         numRegionSegments = polylines.numRegionSegments(i);
%         regionSegments = polylines.regionSegments(i,1:numRegionSegments);
%         CM = sum(X(T(regionSegments,1),:),1)/numRegionSegments;
%         regionPoint = X(T(regionSegments(1),1),:)*0.95 + CM*0.05 ;
%         fprintf(stream,'%d %f %f %f %f %f\n',[regionPoint i sizeRegion]  ); 
%     end
    
    
    fprintf(stream,'\n');
    fclose(stream);
end

function [mesh] = readTriangleOutput( fileName, polylines ,readEdgeMarks, coarseMarkElements)

    checkTriangleOrientation = false;

    nodeName = [fileName '.1.node'];
    elemName = [fileName '.1.ele'];
    edgeName = [fileName '.1.edge'];

    %% read nodes
    fid = fopen(nodeName, 'r');
    outLine = fscanf(fid,'%d %d %d %d\n',4);
    numNodes = outLine(1);
    dim = outLine(2);
    X = fscanf(fid, '%d %f %f %d\n', [(dim+2) numNodes])';
    X = X(:,2:(dim+1));
    fclose(fid);

    %% read elements
    fid = fopen(elemName, 'r');
    outLine = fscanf(fid,'%d %d %d\n',3);
    numElements = outLine(1);
    numNodesElement = outLine(2);
    numAttributesPerElement = outLine(3);
    if(numAttributesPerElement==1)
        elemData = fscanf(fid, '\n %d %d %d %d', [(numNodesElement+2) numElements])';
        elementField = zeros(size(elemData,1),numAttributesPerElement);
        elementField(:,:) = elemData(:,end);
        mesh.elementField = elementField;
    else
        elemData = fscanf(fid, '\n %d %d %d', [(numNodesElement+1) numElements])';
        mesh.elementField = 0;
    end
    T = elemData(:,2:(numNodesElement+1));
    fclose(fid);
    
    %% read edges
    if(readEdgeMarks)
    fid = fopen(edgeName, 'r');
    outLine = fscanf(fid,'%d %d\n',2);
    numEdges = outLine(1);
    numMarks = outLine(2);
    edgeData = fscanf(fid, '%d %d %d %d \n', [4 numEdges])';
    edges = zeros(numEdges,2);
    edges(:,:) = edgeData(:,2:3);
    isBoundaryEdge = edgeData(:,end)>0;
    boundaryMark = edgeData(:,end);
    fclose(fid);
    
    numNodes = size(X,1);
    kTot=numEdges*2;
    i=zeros(kTot,1);
    j=zeros(kTot,1);
    sel=zeros(kTot,1);
    sid=zeros(kTot,1);
    k=1;
    for iedge =1:numEdges
        n1 = edges(iedge,1);
        n2 = edges(iedge,2);
        kn = k+1;
        i(k:kn)=[n1 n2];
        j(k:kn)=[n2 n1];
        if(isBoundaryEdge(iedge))
            %sel(k:kn)=1;
            if(boundaryMark(iedge)==polylines.buildingMark)
                sel(k:kn)=polylines.buildingMark;
            elseif(boundaryMark(iedge)==polylines.boundaryMark)
                sel(k:kn)=polylines.boundaryMark;
            else
                error('not possible')
            end
            
            sid(k:kn) = iedge;
            
%             figure(111)
%             axis equal
%             hold on
%             plot(X(edges(iedge,:),1),X(edges(iedge,:),2),'r')
        else
            sel(k:kn)=-1;
            
            sid(k:kn) = 0;
            
%             figure(111)
%             hold on
%             plot(X(edges(iedge,:),1),X(edges(iedge,:),2),'--k')
        end
        k=kn+1;
    end
%     i=i(1:kn);
%     j=j(1:kn);
%     sel=sel(1:kn);
%     sid=sid(1:kn);
    NNglobal=sparse(i,j,sel,numNodes,numNodes);
    NNedgeId=sparse(i,j,sid,numNodes,numNodes);

    edgeInfo.T = edges(:,1:2);
    edgeInfo.isBoundary = isBoundaryEdge;
    edgeInfo.NN = NNglobal;
    edgeInfo.mark = boundaryMark;
    edgeInfo.buildingMark = polylines.buildingMark;
    edgeInfo.boundaryMark = polylines.boundaryMark;
    edgeInfo.NNedgeId = NNedgeId;
    
    boundaryNodes = unique(edgeInfo.T(find(boundaryMark==polylines.boundaryMark),:));
    mesh.boundaryNodes = boundaryNodes;
    
    else
        edgeInfo = [];
    end
    
    %% Elemental marks for mesh sizing
    elementField = mesh.elementField;
    markedElements = find(elementField);
    coarseElements = markedElements(find(elementField(markedElements)==coarseMarkElements));
    elementField(coarseElements)=1;%elementField(markedElements)=1;
    %sizeMark = ones(size(elementField));
    %sizeMark(coarseElements) = 2;
    %mesh.markSurfaceSize = sizeMark;
    mesh.coarseElements = coarseElements;
    mesh.elementField = elementField;
    
    mesh.numFields = max(elementField);
    
    %% mesh struct
    mesh.T = T;
    mesh.X = X;
    mesh.edges = edgeInfo;
        
    if(checkTriangleOrientation)
        for ielem=1:size(T,1)
            Xelem = X(T(ielem,:),:);
            v1 = [ Xelem(2,:)-Xelem(1,:) 0];
            v2 = [ Xelem(3,:)-Xelem(1,:) 0];
            n = cross(v1,v2);
            if(n(3)<0) 
                disp('Hey! inverted triangle from triangle meshing')
                T(ielem,:) = T(ielem,[1 3 2]);
            end
        end
    end
    
end



