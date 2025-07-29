function [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=getMatrixAdjacentElement_general(...
    T,numNodes,element)

    switch element.type
        
        case{'tri','quad'}
            [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=...
                getMatrixAdjacentElement_2D(T,numNodes,element);
%         case{'tet','hex'}
%             [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=...
%                 getMatrixAdjacentElement_3D(T,numNodes,element);
        case{'tet'}
            [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=...
                getMatrixAdjacentElement_tet(T,numNodes);
            
        case{'hex'}
            [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=...
                getMatrixAdjacentElement_hex(T,numNodes);
        
        case{'pri'}
            [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=...
                getMatrixAdjacentElement_pri(T,numNodes);
            
    end

end

function [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=getMatrixAdjacentElement_2D(...
    T,numNodes,element)
numElements=size(T,1);
numVerticesElement=element.numVertices;
%% Create the global matrix for the search of adjacency
kTot=numElements*numVerticesElement;
i=zeros(kTot,1);
j=zeros(kTot,1);
s=zeros(kTot,1);
k=1;
for iElem=1:numElements
    kn=k+numVerticesElement-1;
    i(k:kn)=iElem;
    j(k:kn)=T(iElem,1:numVerticesElement);
    s(k:kn)=1;
    k=kn+1;
end
ENglobal=sparse(i,j,s,numElements,numNodes);
%%
matrixAdjacentElement=zeros(numElements,element.numEdges);
matrixLocalFaceAdjacentElement=zeros(numElements,element.numEdges);

for iElem=1:numElements
    for iEdge=1:numVerticesElement
        
        nodesEdge = getEdge(T(iElem,:),element,iEdge);
        firstNodeEdge = nodesEdge(1);
        lastNodeEdge  = nodesEdge(length(nodesEdge));
        
        listElements=find(ENglobal(:,firstNodeEdge));       
        theElementsWithTheEdge=listElements(find(ENglobal(listElements,lastNodeEdge)));
%         ENlocal = ENglobal(listElements,:);
%         theElementsWithTheEdge=listElements(find(ENlocal(:,lastNodeEdge)));
%         theAdjacentElement=theElementsWithTheEdge(find(theElementsWithTheEdge~=iElem));
        theAdjacentElement = setdiff(theElementsWithTheEdge,iElem);
        if(isempty(theAdjacentElement))
            matrixAdjacentElement(iElem,iEdge)=0;
            matrixLocalFaceAdjacentElement(iElem,iEdge)=0;
        else
%            T( theAdjacentElement,:)
% iElem
% theAdjacentElement
% T(theAdjacentElement,:)
% T(iElem,:)
            matrixAdjacentElement(iElem,iEdge)=theAdjacentElement;
            matrixLocalFaceAdjacentElement(iElem,iEdge)=...
                find(T(theAdjacentElement,1:numVerticesElement)==lastNodeEdge);
        end
        
    end    
end

end

% function [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=getMatrixAdjacentElement_3D(...
%     T,numNodes,element)
% numElements=size(T,1);
% numVerticesElement=element.numVertices;
% numFacesElement = element.numFaces;
% %% Create the global matrix for the search of adjacency
% kTot=numElements*numVerticesElement;
% i=zeros(kTot,1);
% j=zeros(kTot,1);
% s=zeros(kTot,1);
% k=1;
% for iElem=1:numElements
%     kn=k+numVerticesElement-1;
%     i(k:kn)=iElem;
%     j(k:kn)=T(iElem,1:numVerticesElement);
%     s(k:kn)=1;
%     k=kn+1;
% end
% ENglobal=sparse(i,j,s,numElements,numNodes);
% %%
% matrixAdjacentElement=zeros(numElements,numFacesElement);
% matrixLocalFaceAdjacentElement=zeros(numElements,numFacesElement);
% 
% for iElem=1:numElements
%     for iFace=1:numFacesElement
%         
%         nodesFace = getFace(T(iElem,:), element.order, iFace);
%         
%         listElements=find(ENglobal(:, nodesFace(1) ));       
%         listElements=listElements(find(ENglobal(listElements, nodesFace(2) )));
%         theElementsWithTheFace=listElements(find(ENglobal(listElements, nodesFace(3) )));
%         theAdjacentElement=theElementsWithTheFace(find(theElementsWithTheFace~=iElem));
%         if(isempty(theAdjacentElement))
%             matrixAdjacentElement(iElem,iFace)=0;
%             matrixLocalFaceAdjacentElement(iElem,iFace)=0;
%         else
%             matrixAdjacentElement(iElem,iFace)=theAdjacentElement;
% %             matrixLocalFaceAdjacentElement(iElem,iFace)=...
% %                 find(T(theAdjacentElement,1:numVerticesElement)==lastNodeEdge);
%         end
%         
%     end    
% end
% 
% end