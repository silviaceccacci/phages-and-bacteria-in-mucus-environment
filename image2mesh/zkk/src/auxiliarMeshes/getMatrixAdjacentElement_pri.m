function [matrixAdjacentElement,matrixLocalFaceAdjacentElement] = getMatrixAdjacentElement_pri(T,numNodes)

numElements=size(T,1);
numVerticesElement=6;
%% Create the global matrix for the search of adjacency
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
[ENglobal] = giveConnectivity_ElementToNode(T, numNodes );
%%
numFaces = 5;
matrixAdjacentElement = zeros(numElements,numFaces);
matrixLocalFaceAdjacentElement = zeros(numElements,numFaces);

faceVertices = [1 3 2 0 
                4 5 6 0
                1 2 5 4  
                2 3 6 5
                3 1 4 6 ];
numVerticesFaces = [3 3 4 4 4];
            
for iElem=1:numElements
    
    Telem = T(iElem,:);
    nodesFace = Telem(faceVertices(:,1:3)); % we only need to check the "triangle topology of a quad"
    
    for iFace = 1:numFaces
        if( matrixAdjacentElement(iElem,iFace) ==0)
            listElements=find(ENglobal(:,nodesFace(iFace,1)));
            theElementsWithTheEdge=listElements(find(ENglobal(listElements,nodesFace(iFace,2))));
            theElementsWithTheFace=theElementsWithTheEdge(...
                find(ENglobal(theElementsWithTheEdge,nodesFace(iFace,3))));
            theAdjacentElement=theElementsWithTheFace(find(theElementsWithTheFace~=iElem));
            if(~isempty(theAdjacentElement))
                matrixAdjacentElement(iElem,iFace)=theAdjacentElement;
                v1 = find(T(theAdjacentElement,1:numVerticesElement)==nodesFace(iFace,1));
                v2 = find(T(theAdjacentElement,1:numVerticesElement)==nodesFace(iFace,2));
                v3 = find(T(theAdjacentElement,1:numVerticesElement)==nodesFace(iFace,3));
                [possibleFaces kk] = find(faceVertices==v1) ;
                [indexPosFaces kk] = find(faceVertices(possibleFaces,:)==v2) ;
                possibleFaces = possibleFaces(indexPosFaces);
                [indexPosFaces kk] = find(faceVertices(possibleFaces,:)==v3) ;
                possibleFaces = possibleFaces(indexPosFaces);
                if(length(possibleFaces)==1)
                    matrixLocalFaceAdjacentElement(iElem,iFace)=possibleFaces;
                else
                    error('cagarelapastoret')
                end
                matrixAdjacentElement(theAdjacentElement,possibleFaces) = iElem;
                matrixLocalFaceAdjacentElement(theAdjacentElement,possibleFaces) = iFace;
            end
        end
    end    
end

end
%%
function [matrixAdjacentElement,matrixLocalFaceAdjacentElement] = getMatrixAdjacentElement_tet_proves(T,numNodes)

    %%
    numElements=size(T,1);
    numVerticesElement=element.numVertices;
    %% 
    if(size(T,2)==4)
        faceVertices = [1 2 3 
                        2 1 4
                        2 3 4
                        1 3 4 ];
    else
        error('not implemented');
    end
    kTot=numVerticesElement*numElements;
    i=zeros(kTot,1);
    j=zeros(kTot,1);
    k=zeros(kTot,1);
    con=zeros(kTot,1);
    sel=zeros(kTot,1);
    sed=zeros(kTot,1);
    k=1;
    for iElem=1:numElements
        kn=k+numVerticesElement-1;
        i(k:kn)=(T(iElem,faceVertices(1,:)));
        j(k:kn)=(T(iElem,faceVertices(2,:)));
        k(k:kn)=(T(iElem,faceVertices(3,:)));
        iaux(k:kn)=sort(T(iElem,faceVertices(1,:)));
        jaux(k:kn)=sort(T(iElem,faceVertices(2,:)));
        kaux(k:kn)=sort(T(iElem,faceVertices(3,:)));
        con(k:kn)=1;
        sel(k:kn)=iElem;
        sel(k:kn)=iElem;
        sed(k:kn)=1:numVerticesElement;
        k=kn+1;
    end
    NNglobalCont=sparse(iaux,jaux,kaux,con,numNodes,numNodes,numNodes);
    NNglobal=sparse(i,j,k,sel,numNodes,numNodes,numNodes);
    NNglobal_edgeId=sparse(i,j,k,sed,numNodes,numNodes,numNodes);
    
    list1 = find(NNglobal);
    list2 = find(NNglobalT);
    theList = intersect(list1,list2);
    
    boundaryId = setdiff(list1,theList);
    [in,jn] = ind2sub(size(NNglobal),boundaryId);

    boundaryNodes = unique([in;jn]);
    %%
    matrixAdjacentElement=zeros(numElements,element.numEdges);
    matrixLocalFaceAdjacentElement=zeros(numElements,element.numEdges);

    elems1 = NNglobal(theList);
    elems2 = NNglobalT(theList);
    edges1 = NNglobal_edgeId(theList);
    edges2 = NNglobalT_edgeId(theList);
	linearindex1 = sub2ind(size(matrixAdjacentElement), elems1, edges1);
	%linearindex2 = sub2ind(size(matrixAdjacentElement), elems2, edges2);
    
    matrixAdjacentElement         (linearindex1) = elems2;
    matrixLocalFaceAdjacentElement(linearindex1) = edges2;
    
end

