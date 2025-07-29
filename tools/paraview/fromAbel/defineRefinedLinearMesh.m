function [meshLin] = defineRefinedLinearMesh(mesh,options)
if(isfield(options,'discontinuous')==false)
    options.discontinuous = 1;
end
if(options.discontinuous)
    [meshLin] = defineRefinedLinearMeshDiscont(mesh,options);
else
    [meshLin] = defineRefinedLinearMeshCont(mesh,options);
end
  
end


function [meshLin] = defineRefinedLinearMeshCont(mesh,options)
    activeCostReduction = options.activeCostReduction; % no reduction with this

    if(isfield(options,'exportOrder'))
        refinedOrder = options.exportOrder;
    else
        factorIncreaseOrder = options.factorIncreaseOrder;
        refinedOrder = mesh.element.order*factorIncreaseOrder;
    end
    [meshes]=giveMeshOfDesiredOrder(mesh.X,mesh.T,mesh.element,refinedOrder);
    meshHO = meshes{2};
%% equispace
    numElements = size(mesh.T,1);
    reference_coord = giveReferencePoints(mesh.element);
    V = Vandermonde_LP(mesh.element.order,reference_coord ,mesh.element);
    [L,U,P] = lu(V');
    elementRefined = mesh.element;
    elementRefined.order = refinedOrder;
    elementRefined.distribution='equispaced';
    elementRefined=defineElement(elementRefined);
    
    refined_coord     = giveReferencePoints(elementRefined);
    p = orthopoly2D_general(refined_coord,mesh.element.order,mesh.element);
    shapeFunctions = U\(L\(P*p));
    
    for iElem = 1:numElements 
        meshHO.X(:,meshHO.T(iElem,:)) = mesh.X(:,mesh.T(iElem,:))*shapeFunctions;
    end
%%

    
    linMesh = giveLinearMesh(meshHO.T,meshHO.element);
    elementp1 = setDefaulElement('2D',mesh.element.type,1);
    
    meshLin.X = meshHO.X;
    meshLin.T = linMesh.T;
    meshLin.element = elementp1;
    
    disp('change this on export quads')
    meshLin.quality = ones(size(meshLin.T,1),1);
    meshLin.CM = ones(3,size(meshLin.T,1));
%% rubish
%     %% define linear mesh
%     nDeg = elementRefined.order;
%     switch element.type
%         case 'tri'
%             TLinRef = subtriangulateHOmesh(1:elementRefined.numNod,nDeg);
%         case 'quad'
%             TLinRef.T = giveSubQuadsHighOrderElement(elementRefined);
%     end
%     numSubElementsHOElement = size(TLinRef.T,1);
%     TLinElements =zeros(numElements*numSubElementsHOElement,element.numVertices);
%     CM =zeros(3,numElements*numSubElementsHOElement);
%     qualitySubelements = zeros(numElements*numSubElementsHOElement ,1);
%     if(isfield(mesh,'qualities'))
%        namesQ = fieldnames(mesh.qualities);
%        for i=1:length(namesQ)
%            qualitiesSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
%        end
%     end
%     if(isfield(mesh,'errorDisc'))
%        namesE = fieldnames(mesh.errorDisc);
%        for i=1:length(namesE)
%            errorDisc.(namesE{i}) = zeros(size(qualitySubelements));
%        end
%     end
% 
%     count = 1;
%     count2 = 0;   
%     countHO=0;
%     countLi=0;
%     for iElem = 1:numElements
%         if(activeCostReduction && checkLinearFace(element,mesh.X(:,mesh.T(iElem,:))))
%             count2=count2+1;
%             countFi = count -1 + element.numVertices;
%             X(1:size(mesh.X,1),count:countFi) = mesh.X(:,mesh.T(iElem,1:element.numVertices));
%             TLinElements(count2,:) = count:countFi;
%             qualitySubelements(count2)=mesh.quality(iElem);
%             CM(1:size(mesh.X,1),count2) = mesh.CM(:,iElem);
%             if(isfield(mesh,'qualities'))
%                for i=1:length(namesQ)
%                    qualitiesSubelements.(namesQ{i})(count2) = mesh.qualities.(namesQ{i})(iElem);
%                end
%             end
%             if(isfield(mesh,'errorDisc'))
%                for i=1:length(namesE)
%                    errorDisc.(namesE{i})(count2) = mesh.errorDisc.(namesE{i})(iElem);
%                end
%             end
%             countLi = countLi+1;
%         else
%             countFi = count -1 + nOfRefiNodes;
%             TnewEl = count:countFi;
%             X(:,count:countFi) = mesh.X(:,mesh.T(iElem,:))*shapeFunctions;
%             
%             count2Fi = count2+numSubElementsHOElement;
%             count2 = count2+1;
%             TLinElements(count2:count2Fi,:) = TnewEl(TLinRef.T);
%             
%             if(isfield(mesh,'quality') )
%                 qualitySubelements(count2:count2Fi)=mesh.quality(iElem);
%             end
%             CM(:,count2:count2Fi) = mesh.CM(:,iElem)*ones(1,numSubElementsHOElement);
%             if(isfield(mesh,'qualities'))
%                for i=1:length(namesQ)
%                    qualitiesSubelements.(namesQ{i})(count2:count2Fi) = mesh.qualities.(namesQ{i})(iElem);
%                end
%             end
%             if(isfield(mesh,'errorDisc'))
%                for i=1:length(namesE)
%                    errorDisc.(namesE{i})(count2:count2Fi) = mesh.errorDisc.(namesE{i})(iElem);
%                end
%             end
%             count2 = count2Fi;
%             countHO = countHO+1;
%         end
%         count = countFi+1;
%     end
%     fprintf('HO exported elems = %d\n',countHO)
%     fprintf('Li exported elems = %d\n',countLi)
%     
%     meshLin.X = X(:,1:countFi);
%     meshLin.T = TLinElements(1:count2,:);
%     if(isfield(mesh,'quality'))
%         meshLin.quality = qualitySubelements(1:count2);
%     end
%     meshLin.CM = CM(:,1:count2);
%     if(isfield(mesh,'qualities'))
%        for i=1:length(namesQ)
%            meshLin.qualities.(namesQ{i}) = qualitiesSubelements.(namesQ{i})(1:count2);
%        end
%     end
%     if(isfield(mesh,'errorDisc'))
%        for i=1:length(namesE)
%            meshLin.errorDisc.(namesE{i}) = errorDisc.(namesE{i})(1:count2);
%        end
%     end

end


function [meshLin] = defineRefinedLinearMeshDiscont(mesh,options)

    if(isfield(options,'exportOrder'))
        refinedOrder = options.exportOrder;
    else
        factorIncreaseOrder = options.factorIncreaseOrder;
        refinedOrder = mesh.element.order*factorIncreaseOrder;
    end
    activeCostReduction = options.activeCostReduction;

    element = mesh.element;
    %% define shape functions for the refined mesh
    numElements = size(mesh.T,1);
    reference_coord = giveReferencePoints(mesh.element);
%     V = Vandermonde_LP(mesh.element.order,reference_coord ,mesh.element);
%     [L,U,P] = lu(V');
    elementRefined = mesh.element;
    elementRefined.order = refinedOrder;
    elementRefined.distribution='equispaced';
    elementRefined=defineElement(elementRefined);
    
    refined_coord     = giveReferencePoints(elementRefined);
    nOfRefiNodes = size(refined_coord,1);
%     p = orthopoly2D_general(refined_coord,mesh.element.order,mesh.element);
%     shapeFunctions = U\(L\(P*p));%shapeFunctions  = zeros(nOfNodes,nOfRefiNodes);
    [shapeFunctions]=getShapeFunctions(mesh.element,reference_coord,refined_coord);
    shapeFunctions = shapeFunctions(:,:,1);
    
    X = zeros(3,nOfRefiNodes*numElements);
    
    %% define linear mesh
    nDeg = elementRefined.order;
    switch element.type
        case 'tri'
            TLinRef = subtriangulateHOmesh(1:elementRefined.numNod,nDeg);
        case 'quad'
            TLinRef.T = giveSubQuadsHighOrderElement(elementRefined);
        case 'tet'
            TLinRef.T = giveSubtetsHighOrderElement(elementRefined);
    end
    numSubElementsHOElement = size(TLinRef.T,1);
    TLinElements =zeros(numElements*numSubElementsHOElement,element.numVertices);
    CM =zeros(3,numElements*numSubElementsHOElement);
    qualitySubelements = zeros(numElements*numSubElementsHOElement ,1);
    if(isfield(mesh,'qualities'))
       namesQ = fieldnames(mesh.qualities);
       for i=1:length(namesQ)
           qualitiesSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
       end
    end
    if(isfield(mesh,'distortions'))
       namesQ = fieldnames(mesh.distortions);
       for i=1:length(namesQ)
           distortionsSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
       end
    end
    if(isfield(mesh,'errorDisc'))
       namesE = fieldnames(mesh.errorDisc);
       for i=1:length(namesE)
           errorDisc.(namesE{i}) = zeros(size(qualitySubelements));
       end
    end

    count = 1;
    count2 = 0;   
    countHO=0;
    countLi=0;
    for iElem = 1:numElements
        if(activeCostReduction && checkLinearFace(element,mesh.X(:,mesh.T(iElem,:))))
            count2=count2+1;
            countFi = count -1 + element.numVertices;
            X(1:size(mesh.X,1),count:countFi) = mesh.X(:,mesh.T(iElem,1:element.numVertices));
            TLinElements(count2,:) = count:countFi;
            qualitySubelements(count2)=mesh.quality(iElem);
            CM(1:size(mesh.X,1),count2) = mesh.CM(:,iElem);
            if(isfield(mesh,'qualities'))
               for i=1:length(namesQ)
                   qualitiesSubelements.(namesQ{i})(count2) = mesh.qualities.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'distortions'))
               for i=1:length(namesQ)
                   distortionsSubelements.(namesQ{i})(count2) = mesh.distortions.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'errorDisc'))
               for i=1:length(namesE)
                   errorDisc.(namesE{i})(count2) = mesh.errorDisc.(namesE{i})(iElem);
               end
            end
            countLi = countLi+1;
        else
            countFi = count -1 + nOfRefiNodes;
            TnewEl = count:countFi;
            X(:,count:countFi) = mesh.X(:,mesh.T(iElem,:))*shapeFunctions;
            
            count2Fi = count2+numSubElementsHOElement;
            count2 = count2+1;
            TLinElements(count2:count2Fi,:) = TnewEl(TLinRef.T);
            
            if(isfield(mesh,'quality') )
                qualitySubelements(count2:count2Fi)=mesh.quality(iElem);
            end
            CM(:,count2:count2Fi) = mesh.CM(:,iElem)*ones(1,numSubElementsHOElement);
            if(isfield(mesh,'qualities'))
               for i=1:length(namesQ)
                   qualitiesSubelements.(namesQ{i})(count2:count2Fi) = mesh.qualities.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'distortions'))
               for i=1:length(namesQ)
                   distortionsSubelements.(namesQ{i})(count2:count2Fi) = mesh.distortions.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'errorDisc'))
               for i=1:length(namesE)
                   errorDisc.(namesE{i})(count2:count2Fi) = mesh.errorDisc.(namesE{i})(iElem);
               end
            end
            count2 = count2Fi;
            countHO = countHO+1;
        end
        count = countFi+1;
    end
    fprintf('Elems exported as curved high-order = %d.\n',countHO)
    fprintf('Elems exported as straight-sided    = %d.\n',countLi)
        
    meshLin.X = X(:,1:countFi);
    meshLin.T = TLinElements(1:count2,:);
    if(isfield(mesh,'quality'))
        meshLin.quality = qualitySubelements(1:count2);
    end
    meshLin.CM = CM(:,1:count2);
    if(isfield(mesh,'qualities'))
       for i=1:length(namesQ)
           meshLin.qualities.(namesQ{i}) = qualitiesSubelements.(namesQ{i})(1:count2);
       end
    end
    if(isfield(mesh,'distortions'))
       for i=1:length(namesQ)
           meshLin.distortions.(namesQ{i}) = distortionsSubelements.(namesQ{i})(1:count2);
       end
    end
    if(isfield(mesh,'errorDisc'))
       for i=1:length(namesE)
           meshLin.errorDisc.(namesE{i}) = errorDisc.(namesE{i})(1:count2);
       end
    end

end

%%
%             for iSubElem = 1:numSubElementsHOElement
%                 count2=count2+1;
%                 TLinElements(count2,:) = TnewEl(TLinRef.T(iSubElem,:));
%                 if(isfield(mesh,'quality') )
%                     qualitySubelements(count2)=mesh.quality(iElem);
%                 end
%                 CM(:,count2) = mesh.CM(:,iElem);
%                 if(isfield(mesh,'qualities'))
%                    for i=1:length(namesQ)
%                        qualitiesSubelements.(namesQ{i})(count2) = mesh.qualities.(namesQ{i})(iElem);
%                    end
%                 end
%                 if(isfield(mesh,'errorDisc'))
%                    for i=1:length(namesE)
%                        errorDisc.(namesE{i})(count2) = mesh.errorDisc.(namesE{i})(iElem);
%                    end
%                 end
%             end



%%
% 
% function [meshLin] = defineRefinedLinearMesh(mesh)
% global factorIncreaseOrder;
% 
%     element = mesh.element;
%     %% define shape functions for the refined mesh
%     numElements = size(mesh.T,1);
%     reference_coord = giveReferencePoints(mesh.element);
%     V = Vandermonde_LP(mesh.element.order,reference_coord ,mesh.element);
%     [L,U,P] = lu(V');
%     elementRefined = mesh.element;
%     elementRefined.order = mesh.element.order*factorIncreaseOrder;
%     elementRefined=defineElement(elementRefined);
% %     elementRefined.distribution = 'equispaced';
%     refined_coord     = giveReferencePoints(elementRefined);
%     nOfRefiNodes = size(refined_coord,1);
%     p = orthopoly2D_general(refined_coord,mesh.element.order,mesh.element);
%     shapeFunctions = U\(L\(P*p));%shapeFunctions  = zeros(nOfNodes,nOfRefiNodes);
%     
%     X = zeros(3,nOfRefiNodes*numElements);
%     
%     %% define linear mesh
%     nDeg = elementRefined.order;
%     switch element.type
%         case 'tri'
%             TLinRef = subtriangulateHOmesh(1:elementRefined.numNod,nDeg);
% %                 1:giveNumNodesElementFromOrder_tri(nDeg),nDeg);
%         case 'quad'
%             TLinRef.T = giveSubQuadsHighOrderElement(elementRefined);
% %                 1:giveNumNodesElementFromOrder_quad(nDeg),nDeg);
%     end
%     numSubElementsHOElement = size(TLinRef.T,1);
%     TLinElements =zeros(numElements*numSubElementsHOElement,element.numVertices);
%     CM =zeros(3,numElements*numSubElementsHOElement);
%     qualitySubelements = zeros(numElements*numSubElementsHOElement ,1);
%     if(isfield(mesh,'qualities'))
%        namesQ = fieldnames(mesh.qualities);
%        for i=1:length(namesQ)
%            qualitiesSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
%        end
%     end
%     if(isfield(mesh,'errorDisc'))
%        namesE = fieldnames(mesh.errorDisc);
%        for i=1:length(namesE)
%            errorDisc.(namesE{i}) = zeros(size(qualitySubelements));
%        end
%     end
% 
%     count = 1;
%     count2 = 0;   
%     for iElem = 1:numElements
%         countFi = count -1 + nOfRefiNodes;
%         TnewEl = count:countFi;
%         X(:,count:countFi) = mesh.X(:,mesh.T(iElem,:))*shapeFunctions;
%         for iSubElem = 1:numSubElementsHOElement
%             count2=count2+1;
%             TLinElements(count2,:) = TnewEl(TLinRef.T(iSubElem,:));
%             qualitySubelements(count2)=mesh.quality(iElem);
%             CM(:,count2) = mesh.CM(:,iElem);
%             if(isfield(mesh,'qualities'))
%                for i=1:length(namesQ)
%                    qualitiesSubelements.(namesQ{i})(count2) = mesh.qualities.(namesQ{i})(iElem);
%                end
%             end
%             if(isfield(mesh,'errorDisc'))
%                for i=1:length(namesE)
%                    errorDisc.(namesE{i})(count2) = mesh.errorDisc.(namesE{i})(iElem);
%                end
%             end
%         end
%         count = countFi+1;
%     end
%     
%     meshLin.X = X(:,1:countFi);
%     meshLin.T = TLinElements;
%     meshLin.quality = qualitySubelements;
%     meshLin.CM = CM;
%     if(isfield(mesh,'qualities'))
%        meshLin.qualities = qualitiesSubelements;
%     end
%     if(isfield(mesh,'errorDisc'))
%        meshLin.errorDisc = errorDisc;
%     end
% end


%%




%             if(isfield(mesh,'qualities'))
%                qualitiesSubelements.IdGeo(count2) =  mesh.qualities.IdGeo(iElem);
%                qualitiesSubelements.IdIni(count2) =  mesh.qualities.IdIni(iElem);
%                qualitiesSubelements.ScJac(count2) =  mesh.qualities.ScJac(iElem);
%             end
%        qualitiesSubelements.IdGeo =  zeros(size(qualitySubelements));
%        qualitiesSubelements.IdIni =  zeros(size(qualitySubelements));
%        qualitiesSubelements.ScJac =  zeros(size(qualitySubelements));