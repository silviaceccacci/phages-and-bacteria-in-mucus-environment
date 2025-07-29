function [heightRegion,heighRoofElement]=projectElementWiseFieldAppoximation(mesh,field)

getElementsOverAverage = true;

dirichletBoundary = false;
minimumAtGround = false;

if(minimumAtGround && dirichletBoundary)
    error('one or the other')
end

%forceIntegration = 'gaussPoints';%'pixels','gaussPoints';
global fieldIntegrationType;
forceIntegration =fieldIntegrationType;

if(nargin<3)
    fakeBuildingHeight = 0.0;
end
%%
groundRegion = mesh.groundRegion;

if(isstruct(field))
    field_ground = field.ground;
    field_roof = field.roof;
else
    field_ground = field;
    field_roof = field;
end

%%
[shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement(1,mesh.element.coord,mesh.element.orderQuadrature,mesh.element);

mesh.shapeFunctions = shapeFunctions(:,:,1);
mesh.gaussWeights = gaussWeights;
mesh.gaussPoints = gaussPoints; % we only need the functions, not the derivatives

numNodes = size(mesh.X,1);
Z = zeros(numNodes,1);

if(minimumAtGround)
    Z(:) = inf; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Classify elements into roofs or ground
X2D = mesh.X(:,1:2);

%interestElements = find(mesh.elementField>=groundRegion);% from now on I'll be assuming that the interestElements are the first ones
%interestMesh = mesh;
%interestMesh.T = mesh.T(interestElements,:);

% groundElements = find(mesh.elementField==groundRegion);% from now on I'll be assuming that the interestElements are the first ones
% groundMesh = mesh;
% groundMesh.T = mesh.T(groundElements,:);

roofElements = find(mesh.elementField>groundRegion);% from now on I'll be assuming that the interestElements are the first ones
roofMesh = mesh;
roofMesh.T = mesh.T(roofElements,:);

switch forceIntegration
    case 'pixels';
        disp('projecting...')
%         [projection] = project(field, interestMesh,'fieldToMesh'); % (ngauss,nelem)
%         [projection_roof] = project(field_roof, roofMesh,'fieldToMesh'); % (ngauss,nelem)
%         [projection_ground] = project(field_ground, groundMesh,'fieldToMesh'); % (ngauss,nelem)
%         disp('end projecting...')
% 
%         Zg = projection.z;
%         SFg = projection.shapeF;
%         elemsProj = projection.elements;
%         shapeFProj = projection.shapeF;
%               
%         shapeF_byElement = projection.shapeF_byElement;
%         numPointsInElem = shapeF_byElement.numPointsInElem;
%         shapeF_inElem = shapeF_byElement.shapeF;
%         points_inElem = shapeF_byElement.points;
        
    case 'gaussPoints'
        %disp('projectMeshToField_elemStruct:project only some elements to each field: roof to lidar, ground to topo')
%         xg_g = mesh.X(groundMesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
%                mesh.X(groundMesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
%                mesh.X(groundMesh.T(:,3),1)*mesh.shapeFunctions(3,:);
%         yg_g = mesh.X(groundMesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
%                mesh.X(groundMesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
%                mesh.X(groundMesh.T(:,3),2)*mesh.shapeFunctions(3,:);
%         [Zg_ground] = project(field_ground, [xg_g(:) yg_g(:)],'pointsToField'); % (ngauss,nelem)
%         Zg_ground = reshape(Zg_ground,size(xg_g));
%         
        xg_r = mesh.X(roofMesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
               mesh.X(roofMesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
               mesh.X(roofMesh.T(:,3),1)*mesh.shapeFunctions(3,:);
        yg_r = mesh.X(roofMesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
               mesh.X(roofMesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
               mesh.X(roofMesh.T(:,3),2)*mesh.shapeFunctions(3,:);
        [Zg_roof] = project(field_roof, [xg_r(:) yg_r(:)],'pointsToField'); % (ngauss,nelem)
        Zg_roof = reshape(Zg_roof,size(xg_r));
         
        wg = gaussWeights/2.0;
  
    otherwise
        error('not implemented integration for the force term')
end

%% Averaged roofs
numRegions = mesh.numFields;
heightRegion = zeros(numRegions,1);
regionArea = zeros(numRegions,1);

% holeBuildings = zeros(numRegions,3);
% numNodesBuildings = zeros(numRegions,1);

heighRoofElement = zeros(length(roofElements),1);
areaRoofElement  = zeros(length(roofElements),1);

nodeRegion = zeros(size(mesh.X,1),1);
for locelem = 1:length(roofElements)
    ielem = roofElements(locelem);
    iregion = mesh.elementField(ielem);
    
    elemNodes = mesh.T(ielem,:);
    switch forceIntegration
        case 'pixels'   
            zaux = Zg_roof( points_inElem{ielem} );
            elemHeight = sum(zaux)/numPointsInElem(ielem);%/2.0;
        case 'gaussPoints'
            elemHeight = Zg_roof(locelem,:)*wg;
    end
    elemArea   = computeArea(mesh.X(elemNodes,1:2));

    heightRegion(iregion) = heightRegion(iregion) + elemHeight*elemArea;
    regionArea(iregion)   = regionArea(iregion) + elemArea;
    
    heighRoofElement(locelem) = elemHeight;
    areaRoofElement(locelem)  = elemArea;

    nodeRegion(elemNodes) = iregion;    
        
    %holeBuildings(iregion,1:2) = holeBuildings(iregion,1:2) + sum(mesh.X(elemNodes,1:2),1);
    %numNodesBuildings(iregion) = length(elemNodes);
end
heightRegion = heightRegion./regionArea;

%%==> NO SE QUE CONY ES AIXO... :-(
%epsHeight = 1;
%holeBuildings(:,3) = heightRegion - epsHeight;
%holeBuildings = bsxfun(@rdivide, holeBuildings,numNodesBuildings);

if(getElementsOverAverage)
    heightRegion_new = zeros(numRegions,1);
    regionArea_new   = zeros(numRegions,1);

    epsHeight = 1;
    epsHeightRegions = heightRegion/2;%we will get the elements over the half/third/whatever
    for locelem = 1:length(roofElements)
        ielem = roofElements(locelem);
        iregion = mesh.elementField(ielem);
        
        elemHeight = heighRoofElement(locelem);
        elemArea   = areaRoofElement(locelem);

        %if(elemHeight>heightRegion(iregion)-epsHeight)
        if(elemHeight>epsHeightRegions(iregion))
            heightRegion_new(iregion) = heightRegion_new(iregion) + elemHeight*elemArea;
            regionArea_new(iregion)   = regionArea_new(iregion) + elemArea;
        end
    end
    heightRegion_new = heightRegion_new./regionArea_new;
    
    if(isempty( find(heightRegion_new<heightRegion) ) == false)
        error('this can not be: new heights must be equal or higher')
    end
    
    differenceHeightRegions = abs(heightRegion-heightRegion_new);
    
    fprintf('      Heigh diference after filtering inner courtyards. Min: %4.2f, Max: %4.2f.\n',...
        min(differenceHeightRegions),max(differenceHeightRegions));
    
    heightRegion = heightRegion_new;
end

% roofNodes = find(nodeRegion);
% Z(roofNodes) = heightRegion(nodeRegion(roofNodes)) + fakeBuildingHeight;

%% output
% X3D = [X2D,  Z];
% 
% mesh_facade3D = mesh;
% mesh_facade3D.X = X3D;
% mesh_facade3D.name = [mesh.name '3D'];
% 
% mesh_facade3D.holeBuildings =[];
% %mesh_facade3D.holeBuildings = holeBuildings;
% 
% % figure(234)
% % listP = 1:500:size(field.points,1);
% % plot3(field.points(listP,1),field.points(listP,2),field.z(listP),'*')
% % hold on
% % plot3(X3D(:,1),X3D(:,2),X3D(:,3),'ro')
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