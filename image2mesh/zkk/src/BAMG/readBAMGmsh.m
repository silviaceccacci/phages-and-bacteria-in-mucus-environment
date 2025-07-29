function [mesh] = readBAMGmsh(fileName)

fileNameDat = [fileName '.msh'];
fid = fopen(fileNameDat, 'r');

kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
dim = fscanf(fid,'%d',1)

kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
kk = fgetl(fid);
numNodes = fscanf(fid,'\n %d',1)
X = [fscanf(fid, '\n %f %f %f', [(dim+1) numNodes])];
X = X(1:dim,:);

kk = fscanf(fid, ' \n %s ',1);
numEdges = fscanf(fid, '%d ',1);
kk = [fscanf(fid, '\n %e %e %e', [3 numEdges])];

elemType = fscanf(fid, '\n %s ',1);
if(strcmp(elemType,'Triangles'))
    numNodesElement = 3;
else
    error('not yet included elemType')
end
numElements = fscanf(fid, '%d ',1);
T = [fscanf(fid, '\n %e %e %e %e', [(numNodesElement+1) numElements])]';
T = T(:,1:numNodesElement);

mesh.T = T;
mesh.X = X;

fclose(fid);

% figure(1)
% element = setDefaulElement('2D','tri',1);
% mesh.element = element;
% mesh.X = [mesh.X; zeros(1,size(mesh.X,2))];
% plotMesh_HOElement(mesh.X,mesh.T,mesh.element);

% mesh.fileName = fileName;
% exportMeshParaview(mesh);


