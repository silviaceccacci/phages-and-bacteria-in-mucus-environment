function [ element ] = defineElement( options,phiPointer, DphiPointer, DDphiPointer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of properties of an element
%   dim              dimension (2,3), 
%   numCoord         number of coordinates (2,3)
%   type             ('tri', 'quad', 'tet', 'hex')
%   space            ('plane', 'srf', '3D')
%   order            (1, p)
%   distribution     ('lineal' for order 1,'fekete' or 'equispaced' for HO)
%   numNod           number of nodes of an element
%   numFaces         number of faces (1 for 2D elements, 4 or 6 for 3D)
%   phiPointer, DphiPointer, DDphiPointer (specially for surfaces)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dim = options.dim;
numCoord = options.numCoord;
type = options.type;
order = options.order;
distribution = options.distribution;
if(isfield(options,'quadrature'))
    element.quadrature = options.quadrature;
else
    element.quadrature = 'gaussLobatto';
end
if(isfield(options,'orderQuadrature'))
    element.orderQuadrature = options.orderQuadrature;
end
if(isfield(options,'shapeFunctions'))
    element.shapeFunctions = options.shapeFunctions;
else
   element.shapeFunctions = 'orthopoly'; 
end
%% dim, numCoord, type, space
element.dim = dim ;         % 2, 3
element.numCoord = numCoord ;    % 2, 3
element.type = type;     % 'tri', 'quad', 'tet', 'hex'
if(element.dim == 2 && element.numCoord == 2)
    element.space = 'plane';
elseif(element.dim == 2 && element.numCoord == 3)
    element.space = 'srf' ;
elseif(element.dim == 3 && element.numCoord == 3)
    element.space = '3D';
else
    error('Incorrect element property definition')
end

%% distribution of the nodes
element.order =  order ; % 1 pel linear, p per ordre p
% if( element.order > 1 )
    element.distribution = distribution; % 'fekete'; o 'equispaced';    
% else
%     element.distribution = 'lineal';
% end

%% numNod, numFaces, numEdges, numVertices, numInnerNod, numBoundNod
element.numNod = giveNumNodesElementFromOrder(order,type);
if(          strcmp(type,'tri')  )
    element.numFaces = 1;
    element.numEdges = 3;
    element.numVertices = 3;
    if(order>3)
        element.numInnerNod = giveNumNodesElementFromOrder(order-3,type);
    elseif(order == 3)
        element.numInnerNod = 1;
    else
        element.numInnerNod = 0;
    end
elseif(    strcmp(type,'quad') )
    element.numFaces = 1;
    element.numEdges = 4;
    element.numVertices = 4;
    if(order>2)
        element.numInnerNod = giveNumNodesElementFromOrder(order-2,type);
    elseif(order == 2)
        element.numInnerNod = 1;
    else
        element.numInnerNod = 0;
    end
elseif(    strcmp(type,'tet'))
    element.numFaces = 4;
    element.numEdges = 6;
    element.numVertices = 4;
    if(order>4)
        element.numInnerNod = giveNumNodesElementFromOrder(order-4,type);
    elseif(order == 4)
        element.numInnerNod = 1;
    else
        element.numInnerNod = 0;
    end
elseif(    strcmp(type,'hex'))
    element.numFaces = 6;
    element.numEdges = 12;
    element.numVertices = 8;
    if(order>3)
        element.numInnerNod = giveNumNodesElementFromOrder(order-3,type);
    elseif(order == 3)
        element.numInnerNod = 1;
    else
        element.numInnerNod = 0;
    end
end
element.numBoundNod = element.numNod - element.numInnerNod;

%% Reference element
% if( element.order > 1 )
%     [reference_coord]=giveReferencePoints(order,strcmp(distribution,'fekete'),element);
%     [ shapeFunctions , gaussWeights ]=...
%         computeShapeFunctionsReferenceElement(order,reference_coord, 0, element );
% %     [ shapeFunctionsSubElements , gaussWeightsSubElements , areaIdealSubtriangles ] ...
% %             = getShapeFunctionsOnSubelements(nDeg,reference_coord,element,Winv);
%     element.reference_coord = reference_coord';
%     element.shapeFunctions  = shapeFunctions;
%     element.gaussWeights    = gaussWeights;
% %     element.shapeFunctionsSubElements   = shapeFunctionsSubElements;
% %     element.gaussWeightsSubElements     = gaussWeightsSubElements;
% end


%% Specially for surfaces: parameterization of the surface
if( strcmp(element.space,'srf') && (nargin > 2 || isfield(options,'phiPointer')) )
%     srfNumber = 0; % analytical phi
%     [phiPointer, DphiPointer, DDphiPointer] = choosePointerSurfDefinition(srfNumber);
    if(isfield(options,'phiPointer'))
        element.phiPointer = options.phiPointer;
        element.DphiPointer = options.DphiPointer;
        element.DDphiPointer = options.DDphiPointer;
    else
        element.phiPointer = phiPointer;
        element.DphiPointer = DphiPointer;
        element.DDphiPointer = DDphiPointer;
    end
end

%%

[element] = setGlobalToIJreordering(element);











