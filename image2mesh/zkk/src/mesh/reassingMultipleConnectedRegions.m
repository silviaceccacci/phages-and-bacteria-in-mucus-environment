function [mesh,buildingToBlock] = reassingMultipleConnectedRegions(mesh) 

numNodes = size(mesh.X,1);

EN = mesh.EN;

count = 0;

fieldElements_old = mesh.fieldElements;


mapOldToNewRegions = 1:mesh.numFields;%zeros(mesh.numFields,1);
buildingsInBlock = cell(mesh.numFields,1);

for inode = 1:numNodes
   
    connectedElems = find(EN(:,inode));
    regions = unique(mesh.elementField(connectedElems));
    roofRegions = find(regions>mesh.groundRegion);
    numRoofRegions = sum(roofRegions);
    roofRegions = regions(roofRegions);
    
    if(numRoofRegions>1)
        count = count+1;
        
        theRegion = roofRegions(1);
        buildingsInBlock{theRegion} = [theRegion];
        for i = 2:length(roofRegions)
            iregion = roofRegions(i);
            regionElements = mesh.fieldElements{iregion};
            mesh.fieldElements{iregion} = [];
            mesh.fieldElements{theRegion} = ...
                [mesh.fieldElements{theRegion}; regionElements(:)];
            mesh.elementField(regionElements) = theRegion;
            
            mapOldToNewRegions(iregion) = theRegion;
            buildingsInBlock{theRegion} = [buildingsInBlock{theRegion}  iregion];
        end
    end
    
end

%%

fprintf('   ...Num of buildings: %d\n',mesh.numFields - mesh.groundRegion)
reorderNewRegions = zeros(mesh.numFields,1);
newRegions = mesh.groundRegion;
reorderNewRegions(mesh.groundRegion) = mesh.groundRegion;
newRegionElements = cell(mesh.numFields,1);
for iregion = (1+mesh.groundRegion):mesh.numFields
    if(~isempty(mesh.fieldElements{iregion}))
        newRegions = newRegions + 1;
        newRegionElements{newRegions} = mesh.fieldElements{iregion};
        reorderNewRegions(iregion) = newRegions;
    end
end
mapOldToNewRegions = reorderNewRegions(mapOldToNewRegions);

mesh.fieldElements = newRegionElements(1:newRegions);
mesh.elementField = reorderNewRegions(mesh.elementField);
%mesh.elementField = mapOldToNewRegions(fieldElements_old);
mesh.numFields = newRegions;


% newRegions = mesh.groundRegion;
% mapOldToNewRegions(mesh.groundRegion) = mesh.groundRegion;
% newRegionElements = cell(mesh.numFields,1);
% buildingsInBlock =  cell(mesh.numFields,1);
% for iregion = (1+mesh.groundRegion):mesh.numFields
%     if(~isempty(mesh.fieldElements{iregion}))
%         newRegions = newRegions + 1;
%         elemsInRegion =  mesh.fieldElements{iregion};
%         newRegionElements{newRegions} = elemsInRegion;
%         
%         buildings = mesh.elementField(elemsInRegion);
%         buildingsInBlock{newRegions} = unique(buildings);
%     end
% end
% mesh.fieldElements = newRegionElements(1:newRegions);
% mesh.elementField = mapOldToNewRegions(mesh.elementField);
% mesh.numFields = newRegions;

%buildingToBlock = []
buildingToBlock.numBlocks = newRegions;
buildingToBlock.map = mapOldToNewRegions;
buildingToBlock.blockElems = mesh.fieldElements(1:newRegions);% newRegionElements(1:newRegions);
buildingToBlock.buildingsInBlock = buildingsInBlock(1:newRegions);


fprintf('   ...Num of blocks: %d\n',mesh.numFields-mesh.groundRegion)


end


% function [mesh,buildingToBlock] = reassingMultipleConnectedRegions(mesh) 
% 
% numNodes = size(mesh.X,1);
% 
% EN = mesh.EN;
% 
% count = 0;
% 
% fieldElements_new = cell(mesh.numFields,1);
% elementField_new = ones(size(mesh.elementField)) * mesh.groundRegion;
% 
% mapOldToNewRegions = -ones(mesh.numFields,1);
% mapOldToNewRegions( mesh.groundRegion ) =  mesh.groundRegion;
% 
% buildingsInBlock = cell(mesh.numFields,1);
% 
% for inode = 1:numNodes
%    
%     connectedElems = find(EN(:,inode));
%     regions = unique(mesh.elementField(connectedElems));
%     roofRegions = find(regions>mesh.groundRegion);
%     numRoofRegions = sum(roofRegions);
%     roofRegions = regions(roofRegions);
%     
%     count = count+1;
%         
%     theRegion = count;
%     for i = 1:length(roofRegions)
%         iregion = roofRegions(i);
%         regionElements = mesh.fieldElements{iregion};
%         fieldElements_new{theRegion} = ...
%             [fieldElements_new{theRegion}; regionElements(:)];
%         elementField_new(regionElements) = theRegion;
% 
%         mapOldToNewRegions(iregion) = theRegion;
%         buildingsInBlock{theRegion} = [buildingsInBlock{theRegion}  iregion];
%     end
%     
% %     if(numRoofRegions>1)
% %         count = count+1;
% % 
% %         theRegion = roofRegions(1);
% %         buildingsInBlock{theRegion} = [theRegion];
% %         for i = 2:length(roofRegions)
% %             iregion = roofRegions(i);
% %             regionElements = mesh.fieldElements{iregion};
% %             mesh.fieldElements{iregion} = [];
% %             mesh.fieldElements{theRegion} = ...
% %                 [mesh.fieldElements{theRegion}; regionElements(:)];
% %             mesh.elementField(regionElements) = theRegion;
% %             
% %             mapOldToNewRegions(iregion) = theRegion;
% %             buildingsInBlock{theRegion} = [buildingsInBlock{theRegion}  iregion];
% %         end
% %     end
%     
% end
% 
% %%
% 
% fprintf('   ...Num of buildings: %d\n',mesh.numFields - mesh.groundRegion)
% reorderNewRegions = zeros(mesh.numFields,1);
% newRegions = mesh.groundRegion;
% reorderNewRegions(mesh.groundRegion) = mesh.groundRegion;
% newRegionElements = cell(mesh.numFields,1);
% for iregion = (1+mesh.groundRegion):mesh.numFields
%     if(~isempty(mesh.fieldElements{iregion}))
%         newRegions = newRegions + 1;
%         newRegionElements{newRegions} = mesh.fieldElements{iregion};
%         reorderNewRegions(iregion) = newRegions;
%     end
% end
% mapOldToNewRegions = reorderNewRegions(mapOldToNewRegions);
% 
% mesh.fieldElements = newRegionElements(1:newRegions);
% mesh.elementField = reorderNewRegions(mesh.elementField);
% %mesh.elementField = mapOldToNewRegions(fieldElements_old);
% mesh.numFields = newRegions;
% 
% 
% % newRegions = mesh.groundRegion;
% % mapOldToNewRegions(mesh.groundRegion) = mesh.groundRegion;
% % newRegionElements = cell(mesh.numFields,1);
% % buildingsInBlock =  cell(mesh.numFields,1);
% % for iregion = (1+mesh.groundRegion):mesh.numFields
% %     if(~isempty(mesh.fieldElements{iregion}))
% %         newRegions = newRegions + 1;
% %         elemsInRegion =  mesh.fieldElements{iregion};
% %         newRegionElements{newRegions} = elemsInRegion;
% %         
% %         buildings = mesh.elementField(elemsInRegion);
% %         buildingsInBlock{newRegions} = unique(buildings);
% %     end
% % end
% % mesh.fieldElements = newRegionElements(1:newRegions);
% % mesh.elementField = mapOldToNewRegions(mesh.elementField);
% % mesh.numFields = newRegions;
% 
% %buildingToBlock = []
% buildingToBlock.numBlocks = newRegions;
% buildingToBlock.map = mapOldToNewRegions;
% buildingToBlock.blockElems = mesh.fieldElements(1:newRegions);% newRegionElements(1:newRegions);
% buildingToBlock.buildingsInBlock = buildingsInBlock(1:newRegions);
% 
% 
% fprintf('   ...Num of blocks: %d\n',mesh.numFields-mesh.groundRegion)
% 
% 
% end

