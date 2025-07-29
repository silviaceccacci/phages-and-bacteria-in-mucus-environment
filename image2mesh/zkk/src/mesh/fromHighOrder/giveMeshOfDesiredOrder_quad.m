function [meshes]=giveMeshOfDesiredOrder_quad(X,T,element,desiredOrders,options)
writeMessages = true;
if(nargin>4)
    if(isfield(options,'messages') && options.messages)
        writeMessages = options.messages;
    end    
end
% fekete_selection = 1; % the input mesh is not in fekete configuration -> 0
%                       % the input mesh is in fekete configuration -> use 1
fekete_selection = strcmp(element.distribution,'fekete');
%% Initalize necessary variables
[p_initial] = element.order;

numOrders = length(desiredOrders)+1;
orders = [p_initial desiredOrders];

meshes = cell(numOrders,1) ;
meshes{1}.order = p_initial ; 
meshes{1}.X = X ;
meshes{1}.T = T ;

space_dimension=size(X,1);
numNodes_initial=size(X,2);
numElements=size(T,1);

%% Allocate all the necessary space
for i = 2:numOrders
    iOrder=orders(i);
    
    meshes{i}.order = iOrder ; 
    meshes{i}.X = zeros(space_dimension,1) ;
    meshes{i}.T = zeros(numElements, giveNumNodesElementFromOrder_quad(iOrder) ) ;
end

%% Creation of the meshes of the desired orders
% Create a matrix where connected nodes have value 1
[NNglobal]=giveConnectivity_NodeToNode(T,numNodes_initial);
% Create matrices of neighbouring elemenets
[matrixAdjacentElement,matrixLocalFaceAdjacentElement]=...
    getMatrixAdjacentElement_quad(T,numNodes_initial);

%% Create the new edge nodes
if(writeMessages) disp('creating edge nodes'); end
pairNodes=[ 1 2; 2 3 ; 3 4 ; 4 1];
NNglobal_toMark=NNglobal;
countNewNodes_orders=ones(numOrders,1);
saveNewNameVertices=zeros(numOrders,numNodes_initial);
for iElem=1:numElements %loop on the elements of the mesh
    
    for iLocalNode=1:element.numEdges %loop on each edge
        
        iNode1=T(iElem,pairNodes(iLocalNode,1));
        iNode2=T(iElem,pairNodes(iLocalNode,2));

        if( NNglobal_toMark(iNode1,iNode2) ) %if the edge has not been analysed yet

            % Construct the new T and X matrices for each desired order
            localEdgeNodes=getEdge(T(iElem,:),element,iLocalNode);    
            
            for iOrder=2:numOrders
                
                desiredOrder=orders(iOrder);
                
                vertices=0;
                if(saveNewNameVertices(iOrder,iNode1)==0) 
                    saveNewNameVertices(iOrder,iNode1)=countNewNodes_orders(iOrder)+desiredOrder-1;
                    meshes{iOrder}.X(:, saveNewNameVertices(iOrder,iNode1)) = X(:,iNode1) ;
                else
                    vertices=vertices+1; 
                end
                if(saveNewNameVertices(iOrder,iNode2)==0) 
                    saveNewNameVertices(iOrder,iNode2)=countNewNodes_orders(iOrder)+desiredOrder-vertices;
                    meshes{iOrder}.X(:, saveNewNameVertices(iOrder,iNode2)) = X(:,iNode2) ;
                else
                    vertices=vertices+1; 
                end
                iniVertice = saveNewNameVertices(iOrder,iNode1);
                endVertice = saveNewNameVertices(iOrder,iNode2);
                
                iLocalEdge = getEdge( 1:giveNumNodesElementFromOrder_quad(desiredOrder) , desiredOrder , iLocalNode ) ;
                
                [physicalPoints]= getPointsHigherOrderCurve...
                    (p_initial,localEdgeNodes,X,fekete_selection,desiredOrder);

                meshes{iOrder}.X(:, countNewNodes_orders(iOrder):(countNewNodes_orders(iOrder)+desiredOrder-2))=...
                    physicalPoints(:,2:desiredOrder);
                
                % Assign the new nodes on the edge of the analysed element
                meshes{iOrder}.T(iElem,pairNodes(iLocalNode,1)) = iniVertice;
           
                meshes{iOrder}.T(iElem,iLocalEdge(2:desiredOrder))=...
                    (countNewNodes_orders(iOrder)):(countNewNodes_orders(iOrder)+desiredOrder-2);
              
                meshes{iOrder}.T(iElem,pairNodes(iLocalNode,2)) = endVertice;

                % Assign the new nodes on the edge of the neighbouring element
                if(matrixAdjacentElement(iElem,iLocalNode))
                    
                    neighElem = matrixAdjacentElement(iElem,iLocalNode);
                    
                    iLocalNodeNeigh = matrixLocalFaceAdjacentElement(iElem,iLocalNode);

                    iLocalEdgeNeigh = getEdge( ...
                        1:giveNumNodesElementFromOrder_quad(desiredOrder) , desiredOrder , iLocalNodeNeigh ) ;
                    
                    meshes{iOrder}.T(neighElem,pairNodes(iLocalNodeNeigh,1))= endVertice;

                    meshes{iOrder}.T(neighElem,iLocalEdgeNeigh(2:desiredOrder))=...
                        (countNewNodes_orders(iOrder)+desiredOrder-2):-1:(countNewNodes_orders(iOrder));

                    meshes{iOrder}.T(neighElem,pairNodes(iLocalNodeNeigh,2))= iniVertice;
                                     
                    if(meshes{1}.T(neighElem,pairNodes(iLocalNodeNeigh,2))~=...
                            meshes{1}.T(iElem,pairNodes(iLocalNode,1)))
                        error('one of these elements has an inverted ordering %d or %d',iElem,neighElem);
                    end
                    
                end
                
                %update new elements to add
                countNewNodes_orders(iOrder)=countNewNodes_orders(iOrder)+desiredOrder-vertices+1;
            end
            
            NNglobal_toMark(iNode1,iNode2)=0;
            NNglobal_toMark(iNode2,iNode1)=0;
            
        end % end of the if that checks if the edge had been analysed yet or not
    end % end loop on each edge of the element    
end % end loop on the elements of the mesh

%% Create the inner nodes of the triangles
if(writeMessages) disp('creating inner nodes'); end
feketeDistribution= fekete_selection;
reference_coord= giveReferencePoints(p_initial,feketeDistribution,element);
V = Vandermonde_LP(p_initial,reference_coord,element);
[L,U,P] = lu(V');
nOfNodes = element.numNod;

for iOrder=2:numOrders
    desiredOrder=orders(iOrder);
    optionsElement = element; optionsElement.order = desiredOrder;
    thisOrderElement = defineElement(optionsElement);
    meshes{iOrder}.element = thisOrderElement;

    %Nodes of the reference element
    elementNew = element; elementNew.order = desiredOrder;
    [ elementNew ] = defineElement(elementNew);
    [evalPoints]=giveReferencePoints(elementNew);

    %Evaluation of the shape functions
    numEvaluatedPoints=size(evalPoints,1);
    x=zeros(numEvaluatedPoints,2);
    x(:,:)=evalPoints;

    shapeFunctions = zeros(nOfNodes,numEvaluatedPoints);
%     p = orthopoly2D_deriv_xieta(x,p_initial,elementNew);
    p = orthopoly2D_quad(x,p_initial,element); % assegurar q el q passo aqui esta be! potser es element (vell)
    shapeFunctions(:,:) = U\(L\(P*p));

    numInnerNodes_order= elementNew.numInnerNod;
    
    for iElem=1:numElements
%         iniNodeAux = element.numVertices*desiredOrder+1;
%         fiNodeAux  = element.numVertices*desiredOrder+numInnerNodes_order;
        iniNodeAux = elementNew.numBoundNod+1;
        fiNodeAux  = elementNew.numBoundNod+numInnerNodes_order;
        meshes{iOrder}.T(iElem,iniNodeAux:fiNodeAux) = ...
                countNewNodes_orders(iOrder):(countNewNodes_orders(iOrder)+numInnerNodes_order-1);
            
        newPoints = X(:,T(iElem,:))*shapeFunctions;
        meshes{iOrder}.X(:,meshes{iOrder}.T(iElem,iniNodeAux:fiNodeAux)) = newPoints(:, iniNodeAux:fiNodeAux);

        countNewNodes_orders(iOrder) = countNewNodes_orders(iOrder)+numInnerNodes_order;       
    end
    
end

%% I should now take out the zeros that I allocated
% for iOrder=1:numOrders
%      meshes{iOrder}.X=meshes{iOrder}.X(:,countNewNodes_orders(iOrder));
% end



