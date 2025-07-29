function [boundNodes,boundElements,boundEdges] = giveBoundaryFromConnectivity_tets(...
    T,totalNumNodes,element)

error('do not need this function: call to the function without _tets and eliminate this one')
%      [matrixAdjacentElement,matrixLocalFaceAdjacentElement] =...
%          getMatrixAdjacentElement_tet(T,totalNumNodes,element);
%      [boundElements,boundEdges] = find(matrixAdjacentElement == 0);
%      boundNodes = zeros(totalNumNodes,1);
%      count = 1;
%      for iElem = 1:length(boundElements)
%          % we just add p of the p+1 of the edge, because this way we do not
%          % have to use unique to remove repeated edges
%          count1 = count + element.order -1;
%          nodesEdge = getEdge(T(boundElements(iElem),:),element,boundEdges(iElem));
%          boundNodes(count:count1) = nodesEdge(1:element.order);
%          count = count1 + 1;
%      end
%      boundNodes = boundNodes(1:(count-1));
     
end