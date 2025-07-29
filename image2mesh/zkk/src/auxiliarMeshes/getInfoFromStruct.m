function [ X, T, nDeg, boundaryNodes, boundaryElements, boundaryEdges ] = getInfoFromStruct( mesh )

X = mesh.X;
T = mesh.T;
nDeg = mesh.nDeg;
boundaryNodes = mesh.boundaryNodes;
boundaryElements = mesh.boundaryElements;
boundaryEdges = mesh.boundaryEdges;


if( size(X,2) < 4 ) % if it has dimension numNodes x spaceDim
    
    X = X';
    T = myPermutation(T, nDeg);
    
end