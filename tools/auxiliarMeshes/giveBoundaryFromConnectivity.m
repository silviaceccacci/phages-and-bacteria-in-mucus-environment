function [boundNodes,boundElements,boundEdges] = giveBoundaryFromConnectivity(...
    T,totalNumNodes,element)

  switch element.type 
      case{'tri','quad'}
          [boundNodes,boundElements,boundEdges] = giveBoundaryFromConnectivity_2D(...
                T,totalNumNodes,element);

      case{'tet','hex'}
        [boundNodes,boundElements,boundEdges] = giveBoundaryFromConnectivity_3D(...
            T,totalNumNodes,element);
        
  end

end

function [boundNodes,boundElements,boundFaces] = giveBoundaryFromConnectivity_2D(...
    T,totalNumNodes,element)

     [matrixAdjacentElement,matrixLocalFaceAdjacentElement] =...
         getMatrixAdjacentElement_general(T,totalNumNodes,element);
     [boundElements,boundFaces] = find(matrixAdjacentElement == 0);
     boundNodes = zeros(totalNumNodes,1);
     count = 1;
     for iElem = 1:length(boundElements)
         % we just add p of the p+1 of the edge, because this way we do not
         % have to use unique to remove repeated edges
         count1 = count + element.order ;
         nodesEdge = getEdge(T(boundElements(iElem),:),element,boundFaces(iElem));
         boundNodes(count:count1) = nodesEdge(1:(element.order+1));
         count = count1 + 1;
     end
     boundNodes = unique(boundNodes(1:(count-1)));
     
end


function [boundNodes,boundElements,boundFaces] = giveBoundaryFromConnectivity_3D(...
    T,totalNumNodes,element)

     [matrixAdjacentElement,matrixLocalFaceAdjacentElement] =...
         getMatrixAdjacentElement_general(T,totalNumNodes,element);
     [boundElements,boundFaces] = find(matrixAdjacentElement == 0);
     boundNodes = zeros(totalNumNodes,1);
     count = 1;
     for iElem = 1:length(boundElements)
         boundFaceNodes = getFace(T(boundElements(iElem),:), element, boundFaces(iElem));
         count1 = count + length(boundFaceNodes) -1;
         boundNodes(count:count1) = boundFaceNodes;
         count = count1 + 1;
     end
     boundNodes = unique(boundNodes(1:(count-1)));
     if(size(boundNodes,1)==1)
         boundNodes = boundNodes';
     end
end