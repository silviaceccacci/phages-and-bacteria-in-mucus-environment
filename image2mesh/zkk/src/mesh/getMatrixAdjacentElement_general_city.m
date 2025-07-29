function [matrixAdjacentElement,matrixLocalFaceAdjacentElement,boundaryNodes,NNglobal]=getMatrixAdjacentElement_general_city(...
    T,numNodes,element)
    
    boundaryNodes = [];
    
    NNglobal = [];

    switch element.type
        
        case{'tri','quad'}
            [matrixAdjacentElement,matrixLocalFaceAdjacentElement,boundaryNodes,NNglobal]=...
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
            
    end

end


function [matrixAdjacentElement,matrixLocalFaceAdjacentElement,boundaryNodes,NNglobal]=...
    getMatrixAdjacentElement_2D(...
    T,numNodes,element)
    %%
    numElements=size(T,1);
    numVerticesElement=element.numVertices;
    %% 
    if(size(T,2)==3)
        firstNodes = [ 1 2 3];
        secondNodes = [2 3 1];
    else
        error('not implemented');
    end
    kTot=numVerticesElement*numElements;
    i=zeros(kTot,1);
    j=zeros(kTot,1);
    sel=zeros(kTot,1);
    sed=zeros(kTot,1);
    k=1;
    for iElem=1:numElements
        kn=k+numVerticesElement-1;
        i(k:kn)=T(iElem,firstNodes);
        j(k:kn)=T(iElem,secondNodes);
        sel(k:kn)=iElem;
        sed(k:kn)=1:numVerticesElement;
        k=kn+1;
    end
%     i=i(1:k);
%     j=j(1:kn);
%     sel=sel(1:kn);
%     sed=sed(1:kn);
    NNglobal=sparse(i,j,sel,numNodes,numNodes);
    NNglobal_edgeId=sparse(i,j,sed,numNodes,numNodes);
    NNglobalT=sparse(j,i,sel,numNodes,numNodes);
    NNglobalT_edgeId=sparse(j,i,sed,numNodes,numNodes);
    
    list1 = find(NNglobal);
    list2 = find(NNglobalT);
    theList = intersect(list1,list2);
    
    boundaryId = setdiff(list1,theList);
    %boundaryId2 = setdiff(list2,theList);
    %setdiff(boundaryId,boundaryId2)
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

%% AQUESTA VERSIO FUNCIONA BE I PASSA DE 30 a 5 sec pel cas senzill
% function [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=getMatrixAdjacentElement_2D(...
%     T,numNodes,element)
%     %%
%     numElements=size(T,1);
%     numVerticesElement=element.numVertices;
%     %% 
%     if(size(T,2)==3)
%         firstNodes = [ 1 2 3];
%         secondNodes = [2 3 1];
%     else
%         error('not implemented');
%     end
%     kTot=numVerticesElement*numElements;
%     i=zeros(kTot,1);
%     j=zeros(kTot,1);
%     sel=zeros(kTot,1);
%     sed=zeros(kTot,1);
%     k=1;
%     for iElem=1:numElements
%         kn=k+numVerticesElement-1;
%         i(k:kn)=T(iElem,firstNodes);
%         j(k:kn)=T(iElem,secondNodes);
%         sel(k:kn)=iElem;
%         sed(k:kn)=1:numVerticesElement;
%         k=kn+1;
%     end
%     kn=kn-1;
%     i=i(1:kn);
%     j=j(1:kn);
%     sel=sel(1:kn);
%     sed=sed(1:kn);
%     NNglobal=sparse(i,j,sel,numNodes,numNodes);
%     NNglobal_edgeId=sparse(i,j,sed,numNodes,numNodes);
% %     NNglobalT=sparse(j,i,sel,numNodes,numNodes);
% %     NNglobalT_edgeId=sparse(j,i,sed,numNodes,numNodes);
%     %%
%     disp('Second part getMatrixAdjacentElement')
%     matrixAdjacentElement=zeros(numElements,element.numEdges);
%     matrixLocalFaceAdjacentElement=zeros(numElements,element.numEdges);
% 
%     for iNode=1:numNodes
%         connectedNodes = find( NNglobal( 1:(iNode-1),iNode ) );
%         %connectedNodes = find( NNglobalT( iNode,1:(iNode-1) ) );
%         connectedNodes = connectedNodes(find( NNglobal( iNode, connectedNodes ) ));
%         %connectedNodes = connectedNodes(find( NNglobalT( connectedNodes, iNode ) ));
% 
%         for jloc=1:length(connectedNodes)
%             jNode = connectedNodes(jloc);
%             
%             el1 = NNglobal       (jNode,iNode);
%             ed1 = NNglobal_edgeId(jNode,iNode);
%             el2 = NNglobal       (iNode,jNode);
%             ed2 = NNglobal_edgeId(iNode,jNode);
% 
%             matrixAdjacentElement         (el1,ed1) = el2;
%             matrixLocalFaceAdjacentElement(el1,ed1) = ed2;
%             matrixAdjacentElement         (el2,ed2) = el1;
%             matrixLocalFaceAdjacentElement(el2,ed2) = ed1;
%         end    
%     end
% 
% end
%% with this new version time has been reduced from 30 to 7 (reduced 4 times)
% function [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=getMatrixAdjacentElement_2D(...
%     T,numNodes,element)
%     numElements=size(T,1);
%     numVerticesElement=element.numVertices;
%     %% 
%     if(size(T,2)==3)
%         firstNodes = [ 1 2 3];
%         secondNodes = [2 3 1];
%     else
%         error('not implemented');
%     end
%     kTot=numVerticesElement*numElements;
%     i=zeros(kTot,1);
%     j=zeros(kTot,1);
%     sel=zeros(kTot,1);
%     sed=zeros(kTot,1);
%     k=1;
%     for iElem=1:numElements
%         kn=k+numVerticesElement-1;
%         i(k:kn)=T(iElem,firstNodes);
%         j(k:kn)=T(iElem,secondNodes);
%         sel(k:kn)=iElem;
%         sed(k:kn)=1:numVerticesElement;
%         k=kn+1;
%     end
%     kn=kn-1;
%     i=i(1:kn);
%     j=j(1:kn);
%     sel=sel(1:kn);
%     sed=sed(1:kn);
%     NNglobal=sparse(i,j,sel,numNodes,numNodes);
%     NNglobal_edgeId=sparse(i,j,sed,numNodes,numNodes);
% %     NNglobalT=sparse(j,i,sel,numNodes,numNodes);
% %     NNglobalT_edgeId=sparse(j,i,sed,numNodes,numNodes);
%     %%
%     disp('Second part getMatrixAdjacentElement')
%     matrixAdjacentElement=zeros(numElements,element.numEdges);
%     matrixLocalFaceAdjacentElement=zeros(numElements,element.numEdges);
% 
%     for iNode=1:numNodes
%         connectedNodes = find( NNglobal( 1:(iNode-1),iNode ) );
%         %connectedNodes = find( NNglobalT( iNode,1:(iNode-1) ) );
%         connectedNodes = connectedNodes(find( NNglobal( iNode, connectedNodes ) ));
%         %connectedNodes = connectedNodes(find( NNglobalT( connectedNodes, iNode ) ));
%         
% % %         elem1 = NNglobal(connectedNodes,iNode);
% % %         edge1 = NNglobal_edgeId(connectedNodes,iNode);
% % %         elem2 = NNglobalT(connectedNodes,iNode);
% % %         edge2 = NNglobalT_edgeId(connectedNodes,iNode);
% % %         
% %         elem1 = NNglobal(iNode,connectedNodes);
% %         edge1 = NNglobal_edgeId(iNode,connectedNodes);
% %         elem2 = NNglobalT(iNode,connectedNodes);
% %         edge2 = NNglobalT_edgeId(iNode,connectedNodes);
% % 
% % %         elem1 = NNglobal(connectedNodes,iNode);
% % %         edge1 = NNglobal_edgeId(connectedNodes,iNode);
% % %         elem2 = NNglobal(iNode,connectedNodes);
% % %         edge2 = NNglobal_edgeId(iNode,connectedNodes);
% % 
% %         linearindex1 = sub2ind(size(matrixAdjacentElement), elem1, edge1);
% %         linearindex2 = sub2ind(size(matrixAdjacentElement), elem2, edge2);
% %         
% %         matrixAdjacentElement         (linearindex1) = elem2;
% %         matrixLocalFaceAdjacentElement(linearindex1) = edge2;
% %         matrixAdjacentElement         (linearindex2) = elem1;
% %         matrixLocalFaceAdjacentElement(linearindex2) = edge1;
% 
%         for jloc=1:length(connectedNodes)
%             jNode = connectedNodes(jloc);
%             
%             el1 = NNglobal       (jNode,iNode);
%             ed1 = NNglobal_edgeId(jNode,iNode);
%             el2 = NNglobal       (iNode,jNode);
%             ed2 = NNglobal_edgeId(iNode,jNode);
% % %             el1 = elem1(jloc);
% % %             el2 = elem2(jloc);
% % %             ed1 = edge1(jloc);
% % %             ed2 = edge2(jloc);
%             %if(~(elem2==0))
%             matrixAdjacentElement         (el1,ed1) = el2;
%             matrixLocalFaceAdjacentElement(el1,ed1) = ed2;
%             matrixAdjacentElement         (el2,ed2) = el1;
%             matrixLocalFaceAdjacentElement(el2,ed2) = ed1;
%             %end
%         end    
%     end
% 
% end
%%
% function [matrixAdjacentElement,matrixLocalFaceAdjacentElement]=getMatrixAdjacentElement_2D(...
%     T,numNodes,element)
% 
%     numElements=size(T,1);
%     numVerticesElement=element.numVertices;
%     %% Create the global matrix for the search of adjacency
%     kTot=numElements*numVerticesElement;
%     i=zeros(kTot,1);
%     j=zeros(kTot,1);
%     s=zeros(kTot,1);
%     k=1;
%     for iElem=1:numElements
%         kn=k+numVerticesElement-1;
%         i(k:kn)=iElem;
%         j(k:kn)=T(iElem,1:numVerticesElement);
%         s(k:kn)=1;
%         k=kn+1;
%     end
%     ENglobal=sparse(i,j,s,numElements,numNodes);
%     %%
%     disp('Second part')
%     matrixAdjacentElement=zeros(numElements,element.numEdges);
%     matrixLocalFaceAdjacentElement=zeros(numElements,element.numEdges);
% 
%     isMarked = false(numElements,element.numEdges);
%     
%     for iElem=1:numElements
%         for iEdge=1:numVerticesElement
%             if(~isMarked(iElem,iEdge))
%                 nodesEdge = getEdge(T(iElem,:),element,iEdge);
%                 firstNodeEdge = nodesEdge(1);
%                 lastNodeEdge  = nodesEdge(length(nodesEdge));
% 
%                 listElements=find(ENglobal(:,firstNodeEdge));       
%                 theElementsWithTheEdge=listElements(find(ENglobal(listElements,lastNodeEdge)));
%                 theAdjacentElement = setdiff(theElementsWithTheEdge,iElem);
% 
%                 if(~isempty(theAdjacentElement))
%                     theAdjacentEdge = find(T(theAdjacentElement,1:numVerticesElement)==lastNodeEdge);
%                     
%                     matrixAdjacentElement(iElem,iEdge)=theAdjacentElement;
%                     matrixLocalFaceAdjacentElement(iElem,iEdge)=theAdjacentEdge;
%                     
%                     matrixAdjacentElement(theAdjacentElement,theAdjacentEdge)=iElem;
%                     matrixLocalFaceAdjacentElement(theAdjacentElement,theAdjacentEdge)=iEdge;
%                     
%                     isMarked(theAdjacentElement,theAdjacentEdge) = true;
%                 end
% 
%                 isMarked(iElem,iEdge) = true;
%             end
%         end    
%     end
% 
% end
