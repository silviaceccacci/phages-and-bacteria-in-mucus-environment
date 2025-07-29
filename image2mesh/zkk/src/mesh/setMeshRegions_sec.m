function [mesh] = setMeshRegions(mesh)

    region = mesh.elementField;

    currentRegion = 1;
        
    remainingElements = find(region==0);
        
    while( ~isempty(remainingElements) )
       
        ielem = remainingElements(1);
                
        currentRegion = currentRegion+1;
                
        checkElements = ielem;
        while( ~isempty(checkElements) ) 
%             length(checkElements)
%             currentRegion
            
            region(checkElements) = currentRegion;
            neighbors = mesh.matrixAdjacentElement(checkElements,:);
            neighbors = neighbors(find(neighbors));
            checkElements = neighbors(find(region(neighbors)==0));
        end
        
        remainingElements = remainingElements(find(region(remainingElements)==0));
        
    end

    %currentRegion = currentRegion-1;
    %region = region-1;
    mesh.groundRegion = 1;
    
    mesh.numFields = currentRegion;
    mesh.elementField = region;
    
end

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

%%
% function [mesh] = setMeshRegions(mesh)
% 
%     region = mesh.elementField;
% 
%     currentRegion = 1;
%     
%     regionBoundaries = 0;%zeros(10000,10000); % numRegions,BoundaryEdges
%     
%     listElements = find(region==0);
%     
%     while( ~isempty(listElements) )
%        
%         ielem = listElements(1);
%                 
%         currentRegion = currentRegion+1;
%         
%         region(ielem) = currentRegion;
%         
%         [neighbors,regionBoundaries] = addNotAssignedNeighbors(...
%             [], mesh.matrixAdjacentElement(ielem,:), region,regionBoundaries,currentRegion);
%         while( ~isempty(neighbors) )
%             newNeighbors = [];
%             for jelem=1:length(neighbors)
%                 theElem = neighbors(jelem);
%                 region(theElem) = currentRegion;
%                 [newNeighbors,regionBoundaries] = addNotAssignedNeighbors(...
%                     newNeighbors,mesh.matrixAdjacentElement(theElem,:),region,regionBoundaries,currentRegion);
%             end
%             neighbors = newNeighbors;
%         end
%         
%         listElements = listElements(find(region(listElements)==0));
%         
%     end
% 
%     mesh.numFields = currentRegion;
%     mesh.elementField = region;
%     
% end
