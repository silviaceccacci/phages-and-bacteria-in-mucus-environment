function [element] = setDefaulElement( space, type, order )

    options.type         = type   ; % 'tri', 'quad', 'tet', 'hex' 
    options.order = order;
    options.distribution = 'fekete'; % 'equispaced', 'lineal'
    options.quadrature   = 'gaussLobatto';%'gauss';%'gaussLobatto'
    options.orderQuadrature = options.order + 1 ;
    
    switch space   
        case {'2D','plane'}
            options.dim          = 2       ; % 2, 3
            options.numCoord     = 2       ; % 2, 3

        case {'surface','srf'}
            options.dim          = 2       ; % 2, 3
            options.numCoord     = 3       ; % 2, 3


        case '3D'
            options.dim          = 3       ; % 2, 3
            options.numCoord     = 3       ; % 2, 3       

    end

    [element] = defineElement(options);


end