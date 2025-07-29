function [nv,ne,nf,numElements,numNodes] = giveNumberTypeNodes(mesh)
% this function computes gives the number of different types of nodes of a
% mesh:
%  face nodes
%  inner edge nodes
%  inner vertices
    
    element = mesh.element;
    numElements = size(mesh.T,1);
    numNodes = size(mesh.X,2);

    vertexNodes = unique(reshape(mesh.T(:,1:element.numVertices),element.numVertices*numElements,1));
    nv = length(setdiff(vertexNodes,mesh.boundaryNodes));

    [matrixAdjacentElement,matrixLocalFaceAdjacentElement,EEglobal]=...
                getMatrixAdjacentElement(mesh.T,numNodes);
    [a b] = find(EEglobal);
    numInnerEdges = length(a)/2;
    ne = numInnerEdges * (element.order-1);

    nf = element.numInnerNod * numElements;
end

