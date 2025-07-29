function [heightRegion_topo,heightRegion_lidar]=computeRegionsHeight(mesh,field)

global fieldIntegrationType;
forceIntegration =fieldIntegrationType;

%% PARAMETERS
avoidPatioStrategy = 'overAverage'; %'none'; 'overAverage'; 'percent';
percent = 0.5; % the number of elements taken into account (per evitar integrar patis)

%%
groundRegion = mesh.groundRegion;

if(isstruct(field))
    field_ground = field.ground;
    field_roof = field.roof;
else
    error('Should have two field sources')
    field_ground = field;
    field_roof = field;
end
%%
[shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement(1,mesh.element.coord,mesh.element.orderQuadrature,mesh.element);

mesh.shapeFunctions = shapeFunctions(:,:,1);
mesh.gaussWeights = gaussWeights;
mesh.gaussPoints = gaussPoints; % we only need the functions, not the derivatives
%%
% interestElements = find(mesh.elementField>=groundRegion);% from now on I'll be assuming that the interestElements are the first ones
% interestMesh = mesh;
% interestMesh.T = mesh.T(interestElements,:);

groundElements = find(mesh.elementField==groundRegion);% from now on I'll be assuming that the interestElements are the first ones
groundMesh = mesh;
groundMesh.T = mesh.T(groundElements,:);

roofElements = find(mesh.elementField>groundRegion);% from now on I'll be assuming that the interestElements are the first ones
roofMesh = mesh;
roofMesh.T = mesh.T(roofElements,:);

switch forceIntegration
    case 'pixels';
        disp('projecting...')
        %[projection] = project(field, interestMesh,'fieldToMesh'); % (ngauss,nelem)
        [projection_roof] = project(field_roof, roofMesh,'fieldToMesh'); % (ngauss,nelem)
        [projection_ground] = project(field_ground, groundMesh,'fieldToMesh'); % (ngauss,nelem)
        disp('end projecting...')

        Zg = projection.z;
        SFg = projection.shapeF;
        elemsProj = projection.elements;
        shapeFProj = projection.shapeF;
              
        shapeF_byElement = projection.shapeF_byElement;
        numPointsInElem = shapeF_byElement.numPointsInElem;
        shapeF_inElem = shapeF_byElement.shapeF;
        points_inElem = shapeF_byElement.points;
        
    case 'gaussPoints'

        xg_r = mesh.X(roofMesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
               mesh.X(roofMesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
               mesh.X(roofMesh.T(:,3),1)*mesh.shapeFunctions(3,:);
        yg_r = mesh.X(roofMesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
               mesh.X(roofMesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
               mesh.X(roofMesh.T(:,3),2)*mesh.shapeFunctions(3,:);
        [Zg_ground] = project(field_ground, [xg_r(:) yg_r(:)],'pointsToField'); % (ngauss,nelem)
        Zg_ground = reshape(Zg_ground,size(xg_r));
        [Zg_roof] = project(field_roof, [xg_r(:) yg_r(:)],'pointsToField'); % (ngauss,nelem)
        Zg_roof = reshape(Zg_roof,size(xg_r));
        
        wg = gaussWeights/2.0;
  
    otherwise
        error('not implemented integration for the force term')
end


%% averaged roofs viewed form the lidar and from the topo

switch avoidPatioStrategy
    case 'none'

        [heightRegion_topo,heightRegion_lidar]=computeRegionsMeanElementHeight(...
            roofElements,mesh,Zg_roof,Zg_ground,wg,forceIntegration);

    case 'overAverage'

        [heightRegion_topo,heightRegion_lidar,regionArea,roofElemHeigh_topo,roofElemHeigh_lidar,roofElemArea]=...
            computeRegionsMeanElementHeight(...
                roofElements,mesh,Zg_roof,Zg_ground,wg,forceIntegration);
            
        %heightRegion_lidar_copy = heightRegion_lidar;    
            
        [heightRegion_topo,heightRegion_lidar,regionArea]=computeOverAverageRegionHeight(...
            mesh,roofElements,heightRegion_topo,heightRegion_lidar,regionArea,roofElemHeigh_topo,roofElemHeigh_lidar,roofElemArea);
        
        %[heightRegion_lidar_copy heightRegion_lidar]'
        
        
    otherwise 
        error('not implemented')
end

end

function [heightRegion_topo,heightRegion_lidar,regionArea,roofElemHeigh_topo,roofElemHeigh_lidar,roofElemArea]=...
    computeRegionsMeanElementHeight(...
            roofElements,mesh,Zg_roof,Zg_ground,wg,forceIntegration)
    numRegions = mesh.numFields;
    heightRegion_lidar = zeros(numRegions,1);
    heightRegion_topo = zeros(numRegions,1);
    regionArea = zeros(numRegions,1);
    
    numRoofElements = length(roofElements);
    roofElemHeigh_topo = zeros(numRoofElements,1);
    roofElemHeigh_lidar = zeros(numRoofElements,1);
    roofElemArea = zeros(numRoofElements,1);
    
    for locelem = 1:numRoofElements
        ielem = roofElements(locelem);
        iregion = mesh.elementField(ielem);

        elemNodes = mesh.T(ielem,:);
        switch forceIntegration
            case 'pixels'   
                zaux = Zg_roof( points_inElem{ielem} );
                elemHeight = sum(zaux)/numPointsInElem(ielem);%/2.0;
            case 'gaussPoints'
                %elemHeight_lidar = Zg_roof(ielem,:)*wg;
                %elemHeight_topo = Zg_ground(ielem,:)*wg;
                elemHeight_lidar = Zg_roof(locelem,:)*wg;
                elemHeight_topo = Zg_ground(locelem,:)*wg;
        end
        elemArea   = computeArea(mesh.X(elemNodes,1:2));

        roofElemHeigh_topo(locelem) = elemHeight_topo;
        roofElemHeigh_lidar(locelem) = elemHeight_lidar;
        roofElemArea(locelem) = elemArea;
        
        heightRegion_topo(iregion)  = heightRegion_topo(iregion)  + roofElemHeigh_topo(locelem)*elemArea;
        heightRegion_lidar(iregion) = heightRegion_lidar(iregion) + roofElemHeigh_lidar(locelem)*elemArea;
        regionArea(iregion)         = regionArea(iregion)         + roofElemArea(locelem);
        
        
    end
    heightRegion_topo = heightRegion_topo./regionArea;
    heightRegion_lidar = heightRegion_lidar./regionArea;
end

function [heightRegion_topo_new,heightRegion_lidar_new,regionArea_new]=computeOverAverageRegionHeight(...
    mesh,roofElements,heightRegion_topo,heightRegion_lidar,regionArea,roofElemHeigh_topo,roofElemHeigh_lidar,roofElemArea)
    %parameters:
    percentHeight = 0.9;

    %method:
    numRoofElements = length(roofElements);
    
    heightRegion_lidar_new = zeros(size(heightRegion_lidar));
    regionArea_new = zeros(size(regionArea));
    heightRegion_topo_new = zeros(size(heightRegion_lidar));
    
    for locelem = 1:numRoofElements
        
        ielem = roofElements(locelem);
        iregion = mesh.elementField(ielem);
        heighRegionWithPatio = heightRegion_lidar(iregion);
        
        heightElem_lidar = roofElemHeigh_lidar(locelem);
        heightElem_topo  = roofElemHeigh_topo(locelem);
               
        if(heightElem_lidar > (heighRegionWithPatio*percentHeight) )
        %if(heightElem > (heighRegionWithPatio-1) )
            heightRegion_lidar_new(iregion) = heightRegion_lidar_new(iregion) + heightElem_lidar*roofElemArea(locelem);
            heightRegion_topo_new(iregion) = heightRegion_topo_new(iregion) + heightElem_topo*roofElemArea(locelem);
            regionArea_new(iregion) = regionArea_new(iregion) + roofElemArea(locelem);
        end
        
    end
    
    heightRegion_lidar_new = heightRegion_lidar_new./regionArea_new;
    heightRegion_topo_new  = heightRegion_topo_new./regionArea_new;
    
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