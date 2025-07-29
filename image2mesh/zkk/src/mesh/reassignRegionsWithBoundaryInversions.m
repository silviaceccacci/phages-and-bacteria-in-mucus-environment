function [mesh]=reassignRegionsWithBoundaryInversions(mesh,fields,heightRegion_topo,heightRegion_lidar,tolHeight)

%tolHeight = 5;
[ZGround]=getSampleHeightGround(mesh,fields);

adjElemMat = mesh.matrixAdjacentElement;
elementRegions = mesh.elementField;

numBuildings = 0 ;
for iregion=(mesh.groundRegion+1):mesh.numFields
    
    regionElements = mesh.fieldElements{iregion};
    maxAdjZ = 0;
    for iadj=1:3
        theAdjElems = adjElemMat(regionElements,1);
        groundAdjElems = theAdjElems(find(elementRegions(theAdjElems)==mesh.groundRegion));
        maxAdjZ = max([maxAdjZ max(ZGround(groundAdjElems))]);
    end
    
    diffHeight = heightRegion_lidar(iregion)-maxAdjZ;
    if( diffHeight < tolHeight)
        
        % remove elements from this region and assign them to ground
        regionElements = mesh.fieldElements{iregion};
        mesh.fieldElements{iregion} = [];
        mesh.fieldElements{mesh.groundRegion} = ...
            [mesh.fieldElements{mesh.groundRegion}; regionElements(:)];
        mesh.elementField(regionElements) = mesh.groundRegion;
        
    else
        
        numBuildings = numBuildings+1;
        
    end
    
end
%%

fprintf('   ...Num of initial buildings: %d\n',mesh.numFields - mesh.groundRegion)

mapOldToNewRegions = zeros(mesh.numFields,1);
newRegions = mesh.groundRegion;
mapOldToNewRegions(mesh.groundRegion) = mesh.groundRegion;
newRegionElements = cell(mesh.numFields,1);
for iregion = (1+mesh.groundRegion):mesh.numFields
    if(~isempty(mesh.fieldElements{iregion}))
        newRegions = newRegions + 1;
        newRegionElements{newRegions} = mesh.fieldElements{iregion};
        mapOldToNewRegions(iregion) = newRegions;
    end
end
mesh.fieldElements = newRegionElements(1:newRegions);
mesh.elementField = mapOldToNewRegions(mesh.elementField);
mesh.numFields = newRegions;

fprintf('   ...Num of final   buildings: %d\n',mesh.numFields-mesh.groundRegion)

end

function [ZGround]=getSampleHeightGround(mesh,field)

global fieldIntegrationType;
forceIntegration =fieldIntegrationType;

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

groundElements = find(mesh.elementField==groundRegion);% from now on I'll be assuming that the interestElements are the first ones
groundMesh = mesh;
groundMesh.T = mesh.T(groundElements,:);

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

        xg = mesh.X(mesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
            mesh.X(mesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
            mesh.X(mesh.T(:,3),1)*mesh.shapeFunctions(3,:);
        yg = mesh.X(mesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
            mesh.X(mesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
            mesh.X(mesh.T(:,3),2)*mesh.shapeFunctions(3,:);

        disp('EIII ESTIC ASSIGNANT ZERO A ALGO Q NO ESTA TORBANTTTTT --> canvia al project')
        [Zg_roof] = project(field_roof, [xg(:) yg(:)],'pointsToField'); % (ngauss,nelem)
        [Zg_ground] = project(field_ground, [xg(:) yg(:)],'pointsToField'); % (ngauss,nelem)

        Zg_roof = reshape(Zg_roof,size(xg));
        Zg_ground = reshape(Zg_ground,size(xg));

        wg = gaussWeights/2.0;
  
    otherwise
        error('not implemented integration for the force term')
end

ZGround = max(Zg_ground,[],2);

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