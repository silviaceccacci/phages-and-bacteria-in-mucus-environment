function [ boundaryElements , boundaryEdges ] = getBoundaryElementsFromBoundaryNodes(T, boundaryNodes, element)
 
    numElements = size(T,1);
    nDeg = giveOrderFromConnectivity(T,element);
    
    markedElements = zeros( numElements , element.numVertices);
    
    boundaryElements = zeros( numElements , 1);
    boundaryEdges = zeros( numElements , element.numEdges);
    contBelems = 0;
    
    [ENglobal]=giveConnectivity_ElementToNode(T, max(max(T)) );

    
    for iNode = 1:length(boundaryNodes)
        
        iBnode = boundaryNodes(iNode);
        
        elemsFound = find(ENglobal(:,iBnode));
        
        for iNumElem = 1:length(elemsFound)

            iBelems = elemsFound(iNumElem);
            
            posElems = find(T(iBelems,:)==iBnode);
            
            [isEdge, bEdge] = getEdgeOfNode( posElems , nDeg );
             
             if(isEdge)
                
                if(isEdge && markedElements(iBelems,bEdge) == 0)
                    
                    if(isempty(find(boundaryElements==iBelems)))
                        contBelems = contBelems + 1;
                        boundaryElements(contBelems) = iBelems ; 
                        boundaryEdges(contBelems,bEdge) = 1;
                        markedElements(iBelems,bEdge) = 1;
                    end
                    
                end
                
            end
            
        end
                
        
    end
    
    boundaryElements = boundaryElements(1:contBelems);
    boundaryEdges = boundaryEdges(1:contBelems,:);
