function [ logic ] = checkLinearElement(element,X,options)
% input: X is the coordinates of an element

    if(nargin<3 || strcmp(options.check,'determinant') ) 
        % by constant determinant (for computation)
        [logic] = checkByDeterminant(element,X,options);
    elseif(strcmp(options.check,'edges'))
        % by edges (for visualization)
        [logic] = checkByEdges(element,X,options);
    end

end



function [logic] = checkByEdges(element,X,options)
    logic=true;
    for iEdge = 1:element.numEdges
        Xedge = X(:,getEdge(1:element.numNod,element,iEdge));
        logic = min([logic checkLinearEdge(Xedge,options)]);
        if(logic==false)
            break;
        end
    end        
end

function [logic] = checkByDeterminant(element,X,options)
    logic = true;
    
    if(isfield(options,'tol'))
        tol = options.tol;
    else
        tol = 1e-2;
    end
    
    switch element.type
        case 'tet'
            elementFace = setDefaulElement('srf','tri',element.order);
            
% %             elementFace.quadrature = element.quadrature;
%             elementFace.quadrature = 'gauss';
%             elementFace.orderQuadrature = element.orderQuadrature;
%             [coordRef]=giveReferencePoints_tri(elementFace);
%             coordElem= [-1 -1; 1 -1; -1 1];
%             [shapeFunctions,gaussWeights,gaussPoints]=...
%                 computeShapeFunctionsReferenceElement(element.order,coordRef,coordElem,elementFace);
                        
%             element.vert = vert;
%             [x,w] = getQuadrature(elementFace,N);

            elementFace.shapeFunctions = 'monomial';
    
            [coordRef]=giveReferencePoints_tri(elementFace);
            [shapeFunctions]=getShapeFunctions(elementFace,coordRef,coordRef);
            
        otherwise
            error('not programmed yet')
    end
    
    for iFace = 1:element.numFaces
       
        Tface = getFace(1:size(X,2), element, iFace);
        Xface = X(:,Tface);
        
        nElem = 1;
        nNodElem = size(Tface,2);
        XX = zeros(nElem,1,size(Xface,1),nNodElem);
        XX(1,1,:,:)=Xface;
        Dphi = computeDphi_HO(XX,shapeFunctions);
        detDphi = determinantNDim(Dphi);
        minDet = min(min(min(detDphi)));
        maxDet = max(max(max(detDphi)));
        
        logicFace = abs( (maxDet-minDet)/maxDet ) < tol;
        
        logic = min(logic,logicFace);
        
    end

end