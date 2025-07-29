function [mesh_facade3D]=projectingMeshToField_elemental(mesh,field,fakeBuildingHeight)

groundRegion = mesh.groundRegion;

if(nargin<3)
    fakeBuildingHeight = 50;
end
%%
[shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement(1,mesh.element.coord,mesh.element.orderQuadrature,mesh.element);
mesh.shapeFunctions = shapeFunctions(:,:,1);
mesh.gaussWeights = gaussWeights;
mesh.gaussPoints = gaussPoints; % we only need the functions, not the derivatives

%%

X2D = mesh.X(:,1:2);

% [projection] = project(field, mesh,'gauss'); % (ngauss,nelem)
% %[Zaux] = project(field, mesh,'nodes'); % (ngauss,nelem)
% Zg = projection.z;
% wg = projection.wg;

% here I am assuming that the interestElements are the first ones
interestElements = find(mesh.elementField>=groundRegion);
interestMesh = mesh;
interestMesh.T = mesh.T(interestElements,:);
[projection] = project(field, interestMesh,'gauss'); % (ngauss,nelem)
%[projection] = project(field, interestMesh,'pointsToField'); % (ngauss,nelem)
%[Zaux] = project(field, mesh,'nodes'); % (ngauss,nelem)
Zg = projection.z;
wg = projection.wg;
elemsProj = projection.elements;
shapeFProj = projection.shapeF;

%%
numNodes = size(mesh.X,1);
Z = zeros(numNodes,1);

%% Least-squares minimization for the ground
projectionType_ground = 'minimization'; %minimization,interpolation

if(strcmp(projectionType_ground,'minimization'))
    %M = zeros(numNodes);
    %error('cnstruir matriu esparsa')
    goundElements = interestElements(find(mesh.elementField(interestElements)==groundRegion));
    im = zeros(length(goundElements)*size(mesh.T,2)^2,1);
    jm = zeros(length(goundElements)*size(mesh.T,2)^2,1);
    km = zeros(length(goundElements)*size(mesh.T,2)^2,1);
    mloc1 = [ 1 2 3 1 2 3 1 2 3];
    mloc2 = [ 1 1 1 2 2 2 3 3 3];
    count = 0 ;
    F = zeros(numNodes,1);
    groundNodes = zeros(numNodes,1);
    for locelem = 1:length(goundElements)
        ielem = goundElements(locelem);
        %iregion = mesh.elementField(ielem);
        %if(iregion==groundRegion)
            elemNodes = mesh.T(ielem,:);
            groundNodes(elemNodes) = 1;
            elemArea   = computeArea(mesh.X(elemNodes,1:2));
            
            Mdiag = (mesh.shapeFunctions.*mesh.shapeFunctions)*mesh.gaussWeights;
            Melem = diag(Mdiag);
            Melem(1,2) = (mesh.shapeFunctions(1,:).*mesh.shapeFunctions(2,:))*mesh.gaussWeights;
            Melem(1,3) = (mesh.shapeFunctions(1,:).*mesh.shapeFunctions(2,:))*mesh.gaussWeights;
            Melem(2,3) = (mesh.shapeFunctions(2,:).*mesh.shapeFunctions(3,:))*mesh.gaussWeights;
            Melem(3,2) = Melem(2,3);
            Melem(3,1) = Melem(1,3);
            Melem(2,1) = Melem(1,2);
            Melem = Melem * elemArea;
            
            Felem = bsxfun(@times,mesh.shapeFunctions,Zg(ielem,:))*mesh.gaussWeights;
            Felem = Felem * elemArea;
            F(elemNodes) = F(elemNodes) + Felem;
                        
            %M(elemNodes,elemNodes) = M(elemNodes,elemNodes) + Melem;
            
            im((count+1):(count+length(elemNodes)^2)) = elemNodes(mloc1);
            jm((count+1):(count+length(elemNodes)^2)) = elemNodes(mloc2);
            km((count+1):(count+length(elemNodes)^2)) = Melem(:);
            count = count + length(elemNodes)^2;
        %end
    end
    groundNodes = find(groundNodes);
    
%    im = groundNodes(im);
%?    jm = groundNodes(jm);
    M = sparse(im,jm,km);
    
    M = M(groundNodes,:);
    M = M(:,groundNodes);
    
    F = F(groundNodes);
    Z(groundNodes) = M\F;
    
elseif(strcmp(projectionType_ground,'interpolation'))
    error('dont use this')
    groundNodes = zeros(numNodes,1);
    nodeNeighElem = zeros(size(mesh.X,1),1);
    %for ielem = 1:size(mesh.T,1)
    for locelem = 1:length(interestElements)
        ielem = interestElements(locelem);
        iregion = mesh.elementField(ielem);
        if(iregion==groundRegion)
            elemNodes = mesh.T(ielem,:);
            elemHeight = Zg(ielem,:)*wg;
            Z(elemNodes) = Z(elemNodes) + elemHeight;
            nodeNeighElem(elemNodes) = nodeNeighElem(elemNodes) + 1;

            groundNodes(elemNodes) = 1;
        end
    end
    groundNodes = find(groundNodes);
    Z(groundNodes) = Z(groundNodes)./nodeNeighElem(groundNodes);
         
end
    
    
%% averaged roofs
numRegions = mesh.numFields;
heightRegion = zeros(numRegions,1);
regionArea = zeros(numRegions,1);

nodeRegion = zeros(size(mesh.X,1),1);
%for ielem = 1:size(mesh.T,1)
for locelem = 1:length(interestElements)
    ielem = interestElements(locelem);
    iregion = mesh.elementField(ielem);
    if(iregion>groundRegion)
        elemNodes = mesh.T(ielem,:);
        elemHeight = Zg(ielem,:)*wg;
        %elemHeight = sum(Zg(ielem,:))/length(Zg(ielem,:));
        elemArea   = computeArea(mesh.X(elemNodes,1:2));

        heightRegion(iregion) = heightRegion(iregion) + elemHeight*elemArea;
        %heightRegion(iregion) = Zg(ielem,end);
        %heightRegion(iregion) = Zg(ielem,:)*wg*elemArea;
        regionArea(iregion)   = regionArea(iregion) + elemArea;
        
        nodeRegion(elemNodes) = iregion;
    %elseif(iregion==groundRegion)
            
    end
end

heightRegion = heightRegion./regionArea;

% tuned result
roofNodes = find(nodeRegion);
Z(roofNodes) = heightRegion(nodeRegion(roofNodes)) + fakeBuildingHeight;

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