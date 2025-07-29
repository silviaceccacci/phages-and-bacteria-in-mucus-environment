function [mesh_facade3D,translationVector,translationVector_saved,options,fields,field_topo]...
    =readCityMeshGeometry_STL(fileName,fileNameSTL,tempFolder,...
    meshSize2D_ground,meshAngle2D,edgeMarks,tolHeight,exportStepsToParaview)
%%
options.meshSize = meshSize2D_ground;
options.meshAngle = meshAngle2D;
options.edgeMarks = edgeMarks;
options.tempFolder = tempFolder;

translationVector=[0 0];
translationVector_saved=translationVector;
%%
fprintf('Reading STL geometry...\n')
fv = stlread([fileNameSTL '.stl']);
X = fv.vertices;
T = fv.faces;
clear fv;
%%
fprintf('Removing duplicated nodes from the STL...\n')
[X,T]=removeRepeatedNodes(X',T');
X = X'; T = T';

[X,T]=removeRepeatedNodes(X',T',11);
X = X'; T = T';

mesh_facade3D.X = X;
mesh_facade3D.T = T;
mesh_facade3D.name = fileName;

[consistency,mesh_facade3D]=checkMeshConsistency_usedNodes(mesh_facade3D);

%%
fprintf('Classifying roofs, facades and ground...\n')
[mesh_facade3D]=classifyElementRegionFromSTL(mesh_facade3D,tolHeight);

%% Invent fields
[field_topo]=generateConstantField(0.0);
fields.topo = field_topo;


%%
if(exportStepsToParaview)
    Xsave = mesh_facade3D.X;
    mesh_facade3D.X = bsxfun(@minus,mesh_facade3D.X,[translationVector 0.0]);
    exportTriMeshToParaview(mesh_facade3D)
end
%%
% figure;
% plot3(X(:,1),X(:,2),X(:,3),'*')
% axis equal

end




function [field]=generateConstantField(H)
    
    hx = 1e10;
    hy = hx;
    x0 = [-hx  -hy];
    nx = 2;
    ny = 2;
    
    field.structured = 1;
    field.z = H*ones(nx,ny);
    field.hx = hx;
    field.hy = hy;
    field.x0 = x0;
    field.nx = nx;
    field.ny = ny;
end
