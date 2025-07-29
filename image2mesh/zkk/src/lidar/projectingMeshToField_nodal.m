function [mesh_facade3D]=projectingMeshToField_nodal(mesh,field,fakeBuildingHeight)

groundRegion = mesh.groundRegion;

if(nargin<3)
    fakeBuildingHeight = 50;
end

X2D = mesh.X(:,1:2);
Z = project(field, X2D,'pointsToField');%'nodes');
%X3D = [X2D,  project(field, X2D)];
%X3D(:,3) = X3D(:,3)+mesh.X(:,3);

% averaged roofs
numRegions = mesh.numFields;
heightRegion = zeros(numRegions,1);
regionArea = zeros(numRegions,1);

nodeRegion = zeros(size(mesh.X,1),1);
for ielem = 1:size(mesh.T,1)
    iregion = mesh.elementField(ielem);
    if(iregion>groundRegion)
        elemNodes = mesh.T(ielem,:);
        elemHeight = sum(Z(elemNodes))/length(elemNodes);  
        elemArea   = computeArea(mesh.X(elemNodes,1:2));

        heightRegion(iregion) = heightRegion(iregion) + elemHeight*elemArea;
        regionArea(iregion)   = regionArea(iregion) + elemArea;
        
        nodeRegion(elemNodes) = iregion;
    end
end

heightRegion = heightRegion./regionArea;

% tuned result
%Z(:) = Z(:)+mesh.X(:,3);
for inode=1:size(mesh.X,1)
    iregion = nodeRegion(inode);
    if(iregion>groundRegion)
        Z(inode) = heightRegion(iregion) + fakeBuildingHeight;
    end
end

disp('Here we should ensure that a street node too close to a facane does not get a height value of the facade due to noise in the image..')

%% output

X3D = [X2D,  Z];

mesh_facade3D = mesh;
mesh_facade3D.X = X3D;
mesh_facade3D.name = [mesh.name '3D'];

end

function [area] = computeArea(X)

if(size(X,1)==3)
    v1 = [ X(2,:)-X(1,:) , 0];
    v2 = [ X(3,:)-X(1,:) , 0];
    area = norm(cross(v1,v2))/2.0;
else
    error('not implemented')   
end

end