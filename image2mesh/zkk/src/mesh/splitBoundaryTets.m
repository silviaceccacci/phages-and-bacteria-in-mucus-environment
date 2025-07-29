function [mesh] = splitBoundaryTets(mesh)

    [mesh,mapOldElemToNew]=midNodeSplit(mesh);
    
    [mesh]=nonBoundFacesSplit(mesh,mapOldElemToNew);

end  


function [mesh,mapOldElemToNew]=midNodeSplit(mesh)
    %% Template
%     template = [1 2 3 5
%                 1 2 5 4
%                 2 3 5 4 
%                 3 1 5 4];
    template = [1 2 3 5
                1 4 2 5
                2 4 3 5 
                1 3 4 5];
%     template = [2 4 3 5
%                 1 3 4 5
%                 1 4 2 5 
%                 1 2 3 5];
    %%
    T = mesh.T;
    X = mesh.X;
    numNodes = size(X,1);
    numElems = size(T,1);
    
    % Opcio 1: KK
%     bnodes = zeros(numNodes,1);
%     bnodes(mesh.boundaryNodes) = 1;
%     numBNodesElem = sum(bnodes(T),2);
%     targetElements = find(numBNodesElem==4);
    % Opcio 2: KK a boundaries puc mirar quins elements apareixen dos cops
    % => wrong: I also have to split elements with nodes in the boundary
    % that have no face in the boundary: YES! this happens!!!
%     onesElems = ones(length(mesh.boundaries.elems),1);
%     elemBoundMark = sparse(mesh.boundaries.elems,onesElems,onesElems,numElems,1);
%     targetBFaceElements = find(elemBoundMark>1);
    % Opcio 3: partir qualsevol tet amb tots els nodes al contorn
    isBoundaryNode = zeros(size(mesh.X,1),1);
    isBoundaryNode(mesh.boundaryNodes) = 1;
    targetElements = find(min(isBoundaryNode(mesh.T),[],2)==1);
    
    numTargetElements = length(targetElements);
    Ttarget = [T(targetElements,:) ((numNodes+1):(numNodes+numTargetElements))'];
   
    fprintf('   Number of elements with all nodes in boundary: %d\n',numTargetElements)
    
    mapOldElemToNew=zeros(numElems,4);
    mapOldElemToNew(targetElements,1) = targetElements;
    mapOldElemToNew(targetElements,2) = numElems + (1:numTargetElements);
    mapOldElemToNew(targetElements,3) = numElems+numTargetElements + (1:numTargetElements);
    mapOldElemToNew(targetElements,4) = numElems+2*numTargetElements + (1:numTargetElements);
    
    Tnew1 = Ttarget(:,template(1,:));
    Tnew2 = Ttarget(:,template(2,:));
    Tnew3 = Ttarget(:,template(3,:));
    Tnew4 = Ttarget(:,template(4,:));
    T(targetElements,:) = Tnew1;
    T = [T; Tnew2; Tnew3; Tnew4];
    
    X = [X; (X(Ttarget(:,1),:) + X(Ttarget(:,2),:) + X(Ttarget(:,3),:) + X(Ttarget(:,4),:))/4.0];
    %% get new boundary elements
    numNodesNew = size(X,1);
    EN1 = giveConnectivity_ElementToNode(Tnew1,numNodesNew);
    EN2 = giveConnectivity_ElementToNode(Tnew2,numNodesNew);
    EN3 = giveConnectivity_ElementToNode(Tnew3,numNodesNew);
    EN4 = giveConnectivity_ElementToNode(Tnew4,numNodesNew);
    boundElements = mesh.boundaries.elems;
    boundFaceNodes = mesh.boundaries.faces;
    for e = 1:length(targetElements)
        theElem = targetElements(e);
        ibe = find(boundElements==theElem); % a bit innefficient
        numBfaces = length(ibe);
%         if(numBfaces<=1)
%             numBfaces
%             %error('not all nodes boundary element');
%         end
        for j= 1:numBfaces
            faceNodes = boundFaceNodes(ibe(j),:);
            if(sum(EN1(e,faceNodes))==3)
                newBelem = theElem;
            elseif(sum(EN2(e,faceNodes))==3)
                newBelem = numElems + e;
            elseif(sum(EN3(e,faceNodes))==3)
                newBelem = numElems + numTargetElements+e;
            elseif(sum(EN4(e,faceNodes))==3)
                newBelem = numElems + numTargetElements*2 + e;
            else
                error('not possible')
            end
            boundElements(ibe(j)) = newBelem;
        end
    end
    %%
    mesh.T = T;
    mesh.X = X;
    mesh.boundaries.elems = boundElements;

end

function [mesh]=nonBoundFacesSplit(mesh,mapOldElemToNew)

    %mesh.faces.elements
    %mesh.faces.nodes
        
    isBoundaryNode = zeros(size(mesh.X,1),1);
    isBoundaryNode(mesh.boundaryNodes) = 1;
    
    belemsKK = find(min(isBoundaryNode(mesh.T'))==1);
    if(~isempty(belemsKK))
        mesh.T(belemsKK(1),:)
        isBoundaryNode(mesh.T(belemsKK(1),:))
        belemsKK(1:10)
        error('some elements still with boundary nodes')
    end
    
    splitFaces  = find(min(isBoundaryNode(mesh.innerFaces.nodes)')==1);
    numSplitFaces = length(splitFaces);
    
    fprintf('   Number of non-boundary faces with all nodes in boundary: %d\n',numSplitFaces)
    
    numNewNodes = numSplitFaces;
    newX =(mesh.X(mesh.innerFaces.nodes(splitFaces,1),:)+...
           mesh.X(mesh.innerFaces.nodes(splitFaces,2),:)+...
           mesh.X(mesh.innerFaces.nodes(splitFaces,3),:))/3.0;
    newNodes = size(mesh.X,1) + (1:numNewNodes);
    
    
    numNewElements = 4*length(splitFaces);  % (2+1) + (2+1)
    T = mesh.T;
    Tnew = [T
            zeros(numNewElements,4)];
        
    inew = size(mesh.T,1);
    for iface =1:numSplitFaces
        theFace = splitFaces(iface);
        theFaceNodes = mesh.innerFaces.nodes(theFace,:);
        
        oldElems = mesh.innerFaces.elements(theFace,:);

        for iadj=1:2
            adjElemOld = oldElems(iadj);
            newElements = mapOldElemToNew(adjElemOld,:);
            if(newElements(1)>0) %if the element has been splitted
                %find which is the current adjElement
                for ie=1:length(newElements)
                    adjElement = newElements(ie);
                    isFaceElement = isempty( setdiff(...
                        theFaceNodes,T(adjElement,:)) );
                    if(isFaceElement)
                        break
                    end
                end
                if(~isFaceElement) 
                    error('Something wrong: face element not found')
                end
            else
                adjElement = adjElemOld;
            end
            % Now we have to split adjElement
            Telem = T(adjElement,:);
            noFaceNode = setdiff(Telem,theFaceNodes);
            elemFaceId = find(Telem==noFaceNode);
            
            faceNodesOutNormal = getFace_tet(Telem,elemFaceId);
            n1 = faceNodesOutNormal(1);
            n2 = faceNodesOutNormal(3);
            n3 = faceNodesOutNormal(2);
            n4 = noFaceNode;
            n5 = newNodes(iface);
            
            %if(~isempty(setdiff([n1 n2 n3], theFaceNodes)))
%                 theFaceNodes
%                 Telem
%                 [n1 n2 n3 n4 n5]
%                 isBoundaryNode(Telem)
            %end
            
            TelemsNew =[ n1 n2 n5 n4
                         n2 n3 n5 n4
                         n3 n1 n5 n4 ];
                     
            Tnew(adjElement,:) = TelemsNew(1,:);
            inew=inew+1;
            Tnew(inew,:) = TelemsNew(2,:);
            inew=inew+1;
            Tnew(inew,:) = TelemsNew(3,:);
        end
        
    end
    
    mesh.X = [mesh.X
              newX];
    mesh.T = Tnew(1:inew,:);
          
    %error('aqui hi ha algun errorrrrrrrrrr')
    
%     warning(['In splitBoundaryTets, should also check the faces that'...
%         ' are not in the boundary but have all nodes on boundaries'])

end





