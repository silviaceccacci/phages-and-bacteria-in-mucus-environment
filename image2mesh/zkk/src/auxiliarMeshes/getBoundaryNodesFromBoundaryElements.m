function [ boundaryNodes ] = getBoundaryNodesFromBoundaryElements(T, boundaryElements, boundaryEdges)

nDeg = giveOrderFromConnectivity(T);

boundaryNodes = [];

vecVertex = [1 2 3 1];

for iNumElem = 1:length(boundaryElements)
  
    iBelems = boundaryElements(iNumElem);
    
    for iEdge = 1:3
        
        if( boundaryEdges(iNumElem,iEdge)  == 1)
 
            boundaryNodes = unique( [ boundaryNodes   getEdge(T(iBelems,:), nDeg, iEdge)     T(iBelems,vecVertex(iEdge))     T(iBelems,vecVertex(iEdge+1))     ] );
            
        end
        
    end
    
end