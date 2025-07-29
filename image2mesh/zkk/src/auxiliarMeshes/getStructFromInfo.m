function [ mesh ] = getStructFromInfo(X, T, nDeg, boundaryNodes, boundaryElements, boundaryEdges)

mesh.X = X;
mesh.T = T;
mesh.nDeg = nDeg;
mesh.boundaryNodes = boundaryNodes;

if(nargin<5)
    [ boundaryElements , boundaryEdges ] = getBoundaryElementsFromBoundaryNodes(T, boundaryNodes);    
end

mesh.boundaryElements = boundaryElements;
mesh.boundaryEdges = boundaryEdges;