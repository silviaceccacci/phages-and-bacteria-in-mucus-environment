function [mesh3D]=generate3DCityMesh(fileName,tempFolder,...
    topHeightOverTerrain,mesh_facade3D,...
    translationVector,translationVector_saved,...
    elementOptions,parameters_tet,options,...
    doOptimizeMesh,doSplitBoundaryTets,...
    exportToParaview3D)

tic; fprintf('Build 3D box...\n')
[surfaceMesh]=buildClosedBoxSurfaceMesh(...
    topHeightOverTerrain,mesh_facade3D,elementOptions,parameters_tet,options);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

cm = [translationVector 0.0];%sum(surfaceMesh.X,1)/size(surfaceMesh.X,1);%min(surfaceMesh.X,1);
surfaceMesh.X = bsxfun(@minus,surfaceMesh.X,cm);

% global cmglob;
% cmglob = cm;
% exportCuts(fileName,field_topo,mesh2D,mesh_facade3D)
% 
% error('para aqui ja tinc el q cal')

if(exportToParaview3D)
    exportTriMeshToParaview(surfaceMesh);
end

save(['./output/' fileName '_surface'],'surfaceMesh','mesh_facade3D','parameters_tet','-v7.3')

tic; fprintf('Generating volume mesh...\n')
mesh3D = generateTetMesh(surfaceMesh,parameters_tet,tempFolder);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

save(['./output/' fileName],'mesh3D','cm','-v7.3')


if(doOptimizeMesh)
    tic; fprintf('Optimizing volume mesh...\n')
    mesh3D = optimizeTetMesh(mesh3D,[]);
    elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);
end

mesh3D_notSplitted = mesh3D;
if(doSplitBoundaryTets)
    tic; fprintf('Splitting boundary tets...\n')
    mesh3D = splitBoundaryTets(mesh3D);
    elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);
    
    if(doOptimizeMesh)
        tic; fprintf('Optimizing volume mesh...\n')
        mesh3D = optimizeTetMesh(mesh3D,[]);
        elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);
    end    
end

minD = 1e10;
maxD = 0;
for ielem =1:size(mesh3D.T,1)
    inodes = mesh3D.T(ielem,:);
    Xelem = mesh3D.X(inodes,:);
    A = [Xelem(2,:)-Xelem(1,:)
         Xelem(3,:)-Xelem(1,:)
         Xelem(4,:)-Xelem(1,:)];
    minD = min([minD det(A)]);
    maxD = max([maxD det(A)]);
end

mesh3D.name = [fileName '_3D'];
if(exportToParaview3D)
    exportTetMeshToParaview(mesh3D);
end

mesh3D.X = bsxfun(@plus,mesh3D.X,cm);

save(['./output/' fileName],'mesh3D','mesh3D_notSplitted','translationVector_saved','-v7.3')




