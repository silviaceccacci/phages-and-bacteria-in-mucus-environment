function [mesh2D,options,domainLimitsIni]=generatePlanarCityMesh(fileName,tempFolder,...
    polylines,fields,translationVector,domainLimits,edgeMarks,...
    elementOptions,edgeLengthGround_coarse,meshSize2D_ground,meshAngle2D,...
    tolFacadePoints,tolCloseRegions,minElementAreaAllowed,tolHeight,...
    matlabPlot,exportStepsToParaview,idealizationLevel)

disp('treure el idealizationLevel fora')
if(strcmp(idealizationLevel,'building'))
    fprintf('Meshing 2D at building level\n')
elseif(strcmp(idealizationLevel,'block'))
    fprintf('Meshing 2D at block level\n')
    fprintf('Collapsing unnecessary nodes on contour regions...\n')
else
    idealizationLevel
    error('Not implemented idealization level')
end

tic; 
doCollapseMesh2D=true;
countIt = 0;
doCloseRegions = false;
outScreen = false; optionsParaview.outScreen = outScreen;
maxNumCloseCollapseLoops = 5;
countLoops = 0;
doCheckArea = false;
options.tempFolder = tempFolder;

domainLimitsIni = domainLimits;
outItPraview = true;
while(doCollapseMesh2D || doCloseRegions || doCheckArea)
    countIt = countIt +1;
    fprintf('   ...It: %d (',countIt);
    
    iniNum = size(polylines.T,1);
    
    [polyBox,polyBox2D,domainLimits]=addBoxTo1DGeometry(    ...
        polylines,fields,translationVector,domainLimits,outScreen,edgeLengthGround_coarse);

    fprintf('meshing, ');
    if(~doCheckArea)
        options.meshSize = 0;
        options.meshAngle = 0;
    else
        options.meshSize = meshSize2D_ground;
        options.meshAngle = meshAngle2D;
    end
    options.edgeMarks = edgeMarks;
    mesh = meshConformingPolylines(polyBox,elementOptions,options);

    mesh.area = computeArea2D(mesh);
    minArea = min(mesh.area);

    fprintf('connectivities, ');
    mesh = addMeshConnectivities(mesh);
    mesh.idealizationLevel = idealizationLevel;

    fprintf('regions, ');
    mesh_withRegions = setMeshRegions(mesh);
      
    if(outItPraview)
        mesh_withRegions.name = [fileName '_' int2str(countIt)];
        exportTriMeshToParaview(mesh_withRegions,optionsParaview);
    end
    
    if(~doCheckArea)
        fprintf('collapsing');
        fprintf(' -coll:%d,close:%d-)\n',doCollapseMesh2D,doCloseRegions)
    else
        fprintf('checking too small elements)\n');
        fprintf('      Minimum elemental area ([desired],[allowed]): %1.2e ([%1.2e],[%1.2e])\n',...
            minArea,meshSize2D_ground,minElementAreaAllowed);
    end 
    
    fprintf('   ... collapse boundary of regions: ');
    if(strcmp(idealizationLevel,'building'))
%         [polylines,doCollapseMesh2D,doCloseRegions,countLoops] =...
%             collapseBoundaryOfRegions_buildings(mesh_withRegions,tolFacadePoints,tolCloseRegions,...
%                 doCollapseMesh2D,doCloseRegions,countLoops,maxNumCloseCollapseLoops,...
%                 doCheckArea,minElementAreaAllowed);
         
        %doCheckArea = true
         
         
        doCloseRegions = false;
        doCollapseMesh2D = false;
%         
%         doCollapseMesh2D = true;
%         %doCloseRegions = false;
%         [polylines,doCollapseMesh2D,doCloseRegions,countLoops] =...
%             collapseBoundaryOfRegions(mesh_withRegions,tolFacadePoints,tolCloseRegions,...
%                 doCollapseMesh2D,doCloseRegions,countLoops,maxNumCloseCollapseLoops,...
%                 doCheckArea,minElementAreaAllowed);
%         doCloseRegions = false;
    elseif(strcmp(idealizationLevel,'block'))
        [polylines,doCollapseMesh2D,doCloseRegions,countLoops] =...
            collapseBoundaryOfRegions(mesh_withRegions,tolFacadePoints,tolCloseRegions,...
                doCollapseMesh2D,doCloseRegions,countLoops,maxNumCloseCollapseLoops,...
                doCheckArea,minElementAreaAllowed);
    end
    
    doCheckArea = ~doCollapseMesh2D && ~doCloseRegions && ~doCheckArea;
        
    %if(maxNumCloseCollapseLoops<countIt)
    %    doCheckArea=false; doCollapseMesh2D=false; doCloseRegions=false;
    %end
    
    fprintf('      Initial/collapsed boundaries: %d/%d\n',iniNum,size(polylines.T,1));
end
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

if(strcmp(idealizationLevel,'building'))
    fprintf('  ...number of blocks: %i\n',mesh_withRegions.numBlocks);
    fprintf('  ...number of buildings: %i\n',mesh_withRegions.numBuildings);
end

tic; fprintf('Adding boundary box...\n')
[polyBox,polyBox2D,domainLimits]=addBoxTo1DGeometry(...
    polylines,fields,translationVector,domainLimits,outScreen,edgeLengthGround_coarse);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

if(domainLimits.sizeBox)
    polyBox.sizeBox = domainLimits.sizeBox;
    polyBox.sizeOut = meshSize2D_ground_coarse;
    polyBox.coordIn = [min(polyBox.X(:,1))+0.1 min(polyBox.X(:,2))+0.1];

    tic; fprintf('Adding coarse boundary box...\n')
    domainLimits.type = domainLimits.typeOut;
    [polyBox,polyBox2D,domainLimits]=addBoxTo1DGeometry(...
        polyBox,fields,translationVector,domainLimits,outScreen,edgeLengthGround_coarse);
    elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);
    polyBox.coordOut = [min(polyBox.X(:,1))+0.1 min(polyBox.X(:,2))+0.1];
end

numFields = size(mesh_withRegions.fieldElements,1);
regionsToMesh.position = zeros(numFields-mesh_withRegions.groundRegion,2);
regionsToMesh.mark = zeros(numFields-mesh_withRegions.groundRegion,1);
for ifield=(mesh_withRegions.groundRegion+1):numFields
    elem = mesh_withRegions.fieldElements{ifield}(1);
    point = mesh_withRegions.X(mesh_withRegions.T(elem,:),:);
    regionsToMesh.position(ifield-1,:) = sum(point,1)/length(mesh_withRegions.T(elem,:));
    regionsToMesh.mark(ifield-1) = ifield;%0;
end

tic; fprintf('Meshing...\n')
options.meshSize = meshSize2D_ground;
options.meshAngle = meshAngle2D;
options.edgeMarks = edgeMarks;
polyBox.extraRegions = regionsToMesh;
mesh = meshConformingPolylines(polyBox,elementOptions,options);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

[areaElements] = computeArea2D(mesh);
minArea = min(areaElements);
fprintf('   Minimum elemental area ([desired],[allowed]): %e ([%e],[%e])\n',minArea,meshSize2D_ground,minElementAreaAllowed);

%     mesh.name='putaMerda';
%     exportTriMeshToParaview(mesh,optionsParaview);
    
tic; fprintf('Setting mesh info (connectivity, neighbors)...\n')
mesh = addMeshConnectivities(mesh);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

mesh.idealizationLevel = idealizationLevel;

tic; fprintf('Finding regions...\n')
mesh_withRegions = setMeshRegions(mesh);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

if(outItPraview)
    mesh_withRegions.name = [fileName '_' int2str(countIt+1)];
    exportTriMeshToParaview(mesh_withRegions,optionsParaview);
end
if(size(polylines.X,1)>0)
    polylines.X = bsxfun(@plus,polylines.X,translationVector);
    polyBox2D.X = bsxfun(@plus,polyBox2D.X,translationVector);
    mesh_withRegions.X = bsxfun(@plus,mesh_withRegions.X,translationVector);
end
  


if(strcmp(idealizationLevel,'block'))
    
    tic; fprintf('Reassign not building regions to ground...\n')
    [heightRegion_topo,heightRegion_lidar]=computeRegionsHeight(mesh_withRegions,fields);
    [mesh_withRegionsNew1]=reassignRegionsThroughHeights(...
        mesh_withRegions,heightRegion_topo,heightRegion_lidar,tolHeight);
    elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

    %[mesh_withRegionsNew2]=reassignRegionsWithBoundaryInversions(...
    %    mesh_withRegionsNew1,fields,heightRegion_topo,heightRegion_lidar,tolHeight);
    mesh_withRegionsNew2 = mesh_withRegionsNew1;

    tic; fprintf('Reassign multiple connected regions...\n')
    % if two regions are connected through a node, and we are working at a "illla" level,
    % they must be the same building (required for topological construction)
    [mesh_withRegionsNew3] = reassingMultipleConnectedRegions(mesh_withRegionsNew2); 
    elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);
    
elseif(strcmp(idealizationLevel,'building'))
    mesh_withRegionsNew3 = mesh_withRegions;
%     tic; fprintf('Reassign multiple connected regions...\n')
%     [meshKK,buildingToBlock] = reassingMultipleConnectedRegions(mesh_withRegionsNew2); 
%     elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);
%     mesh_withRegionsNew3 = mesh_withRegionsNew2;
%     mesh_withRegionsNew3.buildingToBlock = buildingToBlock;
end

mesh_withRegionsNew3.name = fileName;
if(exportStepsToParaview)
    Xsave = mesh_withRegionsNew3.X;
    mesh_withRegionsNew3.X = bsxfun(@minus,mesh_withRegionsNew3.X,translationVector);
    exportTriMeshToParaview(mesh_withRegionsNew3)
    mesh_withRegionsNew3.X = Xsave;
end

mesh2D = mesh_withRegionsNew3;

if(matlabPlot)
    fprintf('Plotting...\n');    figure(3);
    plotMesh(mesh2D.X,mesh2D.T,mesh2D.elementField)
    caxis([min(mesh2D.elementField) max(mesh2D.elementField)])
end
