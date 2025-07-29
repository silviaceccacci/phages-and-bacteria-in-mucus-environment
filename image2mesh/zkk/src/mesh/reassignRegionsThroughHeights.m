function [mesh]=reassignRegionsThroughHeights(mesh,heightRegion_topo,heightRegion_lidar,tolHeight)

% remeshingRegions.regions = cell(mesh.numFields,1);
% remeshingRegions.numRegions = 0;

numBuildings = 0 ;
for iregion=(mesh.groundRegion+1):mesh.numFields
   
    diffHeight = heightRegion_lidar(iregion)-heightRegion_topo(iregion);
    if( diffHeight < tolHeight)
        
        % remove elements from this region and assign them to ground
        regionElements = mesh.fieldElements{iregion};
        regionNodes = unique(mesh.T(regionElements,:));
        mesh.fieldElements{iregion} = [];
        mesh.fieldElements{mesh.groundRegion} = ...
            [mesh.fieldElements{mesh.groundRegion}; regionElements(:)];
        mesh.elementField(regionElements) = mesh.groundRegion;
        
%         remeshingRegions.numRegions = remeshingRegions.numRegions+1;
%         regionStruct.xcm = mesh.X(regionNodes,1)/length(regionNodes);
%         regionStruct.ycm = mesh.X(regionNodes,2)/length(regionNodes);
%         remeshingRegions.regions{remeshingRegions.numRegions} = regionStruct;
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