function [mesh] = setMeshRegions(mesh)

    mesh.groundRegion = 1; % to change this change triangle meshing marks
        
    if(~isfield(mesh,'numFields') || mesh.numFields==1)
        region = mesh.elementField;

        currentRegion = mesh.groundRegion;

        remainingElements = find(region==0);

        groundElements = find(region==mesh.groundRegion);
        regionElements = cell(length(remainingElements),1);
        regionElements{mesh.groundRegion} = groundElements;

        while( ~isempty(remainingElements) )

            ielem = remainingElements(1);

            currentRegion = currentRegion+1;

            checkElements = ielem;
            while( ~isempty(checkElements) ) 
                region(checkElements) = currentRegion;
                regionElements{currentRegion} = [regionElements{currentRegion}; checkElements(:)];

                neighbors = mesh.matrixAdjacentElement(checkElements,:);
                neighbors = neighbors(find(neighbors));
                checkElements = unique(neighbors(find(region(neighbors)==0)));
            end

            remainingElements = remainingElements(find(region(remainingElements)==0));

        end

        mesh.numFields = currentRegion;
        mesh.elementField = region;
        mesh.fieldElements = regionElements(1:currentRegion); 

        if(isfield(mesh,'idealizationLevel') && strcmp(mesh.idealizationLevel,'building'))    
            
            numBlocks = max(mesh.elementField)-min(mesh.elementField);
            mesh.elementField_blocks = mesh.elementField;
            
            [mesh] = setMeshSubregions(mesh);
            
            numBuildings = max(mesh.elementField)-min(mesh.elementField);
            
            mesh.numBlocks = numBlocks;
            mesh.numBuildings = numBuildings;
%             fprintf('      ...number of blocks: %i\n',numBlocks);
%             fprintf('      ...number of buildings: %i\n',numBuildings);
        end

    else
        numFields = mesh.numFields;
        regionElements = cell(numFields,1);
        for ifield = 1:numFields
            regionElements{ifield} = find(mesh.elementField==ifield);
        end        
        mesh.fieldElements = regionElements; 
    end
end


function [mesh] = setMeshSubregions(mesh)

    fprintf(' setting subregions, ');
    edgeInfo = mesh.edges;
    NN = edgeInfo.NN;

    currentRegion = mesh.groundRegion;
    region = mesh.elementField;
    
    %newRegion = region;
    newRegion = ones(size(region)) * inf;
    groundElements = find(region==mesh.groundRegion);
    newRegion(groundElements) = mesh.groundRegion;
    
    isMarked = false(size(mesh.T,1),1);
    isMarked(groundElements) = true;
    
    elemEdges = [1 2; 2 3 ; 3 1];
    
    buildingsInBlock = cell(mesh.numFields,1);
    
    for iregion = (mesh.groundRegion+1):mesh.numFields
        regionElements = mesh.fieldElements{iregion};
        
        remainingElements = regionElements;

        
        while( ~isempty(remainingElements) )
            ielem = remainingElements(1);
            currentRegion = currentRegion+1;
            buildingsInBlock{iregion} = [buildingsInBlock{iregion} currentRegion];

            checkElements = ielem;
            while( ~isempty(checkElements) ) 
                
                newRegion(checkElements) = currentRegion;
                isMarked(checkElements) = true;
                
                remainingElements = setdiff(remainingElements,checkElements);
                
                neighbors = mesh.matrixAdjacentElement(checkElements,:);
           
                newCheckElements = [];
                for jelem = 1:length(checkElements)
                    currentElement = checkElements(jelem);
                    
                    for ineigh = 1:size(neighbors,2)%length(neighbors)

                        if(neighbors(jelem,ineigh)~=0)
                            neighElem = neighbors(jelem,ineigh);
                            
                            edgeNodes = mesh.T(currentElement,elemEdges(ineigh,:));
                            if( NN(edgeNodes(1),edgeNodes(2)) == -1 ) % is not on geometry edge
                                if(~isMarked(neighElem))
                                    newCheckElements = [newCheckElements neighElem];
                                end
                            end
                        end
                    end
                end
                checkElements = unique(newCheckElements);
            end
        end
        
    end
        
    buildingToBlock.numBlocks = mesh.numFields;
    [buildUnique, index] = unique(newRegion);
    mapBuildingToBlock = mesh.elementField(index);
    buildingToBlock.map = mapBuildingToBlock;
    buildingToBlock.blockElems = mesh.fieldElements;
    buildingToBlock.buildingsInBlock = buildingsInBlock;
    mesh.buildingToBlock = buildingToBlock;
    
    mesh.numFields = currentRegion;
    mesh.elementField = newRegion;
    %mesh.fieldElements =[]; % REDO!!!!!
    
    numFields = mesh.numFields;
    regionElements = cell(numFields,1);
    for ifield = 1:numFields
        regionElements{ifield} = find(mesh.elementField==ifield);
    end        
    mesh.fieldElements = regionElements; 
    
end



%%
function [neighbors,regionBoundaries]=addNotAssignedNeighbors(neighbors,neighElem,region,regionBoundaries,currentRegion)
    
    neighElem = neighElem(find(neighElem));
    neighbors = [neighbors neighElem(find(region(neighElem)==0))];
    
%     for ielem = 1:length(neighElem)
%        theElem = neighElem(ielem);
%        if(theElem>0)
%            if(region(theElem)==0)
%                neighbors = [neighbors theElem];
%            %elseif(region(theElem~=currentRegion)
%            %    
%            end
%        end
%     end

end

function [mesh] = assignRegionElements(mesh)
error('TOO SLOW: included directly in setMeshRegions')
    fprintf('   ...building data structure:');
    percentage = 10;
    numMod = floor(size(mesh.T,1)/percentage);
    acumPercentage = 0;

    regionElements = cell(mesh.numFields,1);
    for ielem = 1:size(mesh.T,1)
        if(mod(ielem,numMod)==0)
            acumPercentage = acumPercentage+percentage;
            fprintf(' %3.0d%%',acumPercentage)
        end
        iregion = mesh.elementField(ielem);
        regionElements{iregion} = [regionElements{iregion} ielem];
    end
    
    mesh.fieldElements = regionElements;
    
    fprintf('\n');
end