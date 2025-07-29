function [mesh] = removeRegionsAsRoofs(mesh,regionsToRemove) 

%numNodes = size(mesh.X,1);

groundRegion = mesh.groundRegion;

%EN = mesh.EN;


if(mesh.specialBuildings.doTreat)
    specialRegions = zeros(length(regionsToRemove),2);
end


for iregion = 1:length(regionsToRemove)
   
    theRegion = regionsToRemove(iregion);
    
    regionElements = mesh.fieldElements{theRegion};
    mesh.fieldElements{theRegion} = [];
    mesh.fieldElements{groundRegion} = ...
        [mesh.fieldElements{groundRegion}; regionElements(:)];
    mesh.elementField(regionElements) = groundRegion;
    
    if(mesh.specialBuildings.doTreat)
        regionNodes = unique(mesh.T(regionElements,:));
        specialRegions(:,iregion) = sum(mesh.X(regionNodes,:),1)/length(regionNodes);
        % aqui hauria de guardar tamb els elements de cada regio especial...
        % els hauria de tenir amb una id diferent...
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
% mesh.fieldElements = newRegionElements(1:newRegions);
% mesh.elementField = mapOldToNewRegions(mesh.elementField);
% mesh.numFields = newRegions;
mesh.fieldElements = newRegionElements(1:newRegions);
%nonFacadeElements = find(mesh.elementField>=groundRegion);
%mesh.elementField = mapOldToNewRegions(mesh.elementField(nonFacadeElements));
mesh.elementField = mapOldToNewRegions(mesh.elementField);
mesh.numFields = newRegions;


fprintf('   ...Num of final   buildings: %d\n',mesh.numFields-mesh.groundRegion)

if(mesh.specialBuildings.doTreat)
    if(isfield(mesh,'specialRegions'))
        specialRegions = [mesh.specialRegions; specialRegions];
    end
    mesh.specialRegions = specialRegions;
end
end