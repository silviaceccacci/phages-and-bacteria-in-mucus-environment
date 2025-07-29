function [mesh]=recomputeGroundMarks(mesh)

     MAE = getMatrixAdjacentElement_general(mesh.T,size(mesh.X,1),mesh.element);

     regions = mesh.elementField;
          
     groundElements = mesh.extendedElementList';
     length(groundElements)
     groundElements = find(regions==mesh.groundRegion);
     length(groundElements)
     
     numGroundElems = length(groundElements);
     numGroundElems_pre = 0;
     while(numGroundElems~=numGroundElems_pre)
         numGroundElems_pre = numGroundElems;
         groundElements_new = groundElements;
         for iedge =1:3
             listNeighbors = MAE(groundElements,iedge); %find neighbors of ground elements
             listNeighbors = listNeighbors(find(listNeighbors)); %find actual neigbors (not boundary)
             listNeighbors = listNeighbors( find( regions(listNeighbors)==mesh.roofRegion ) ) ; %find neighors that are roofs (not possible)
             groundElements_new = unique([groundElements_new; listNeighbors]); %roof neighbors are aactually ground
         end
         groundElements = groundElements_new;
         numGroundElems = length(groundElements);
         regions(groundElements) = mesh.groundRegion;
         
%          [numGroundElems_pre numGroundElems]
     end
     
     mesh.elementField = regions;

     mesh = rmfield(mesh,'extendedElementList');
end