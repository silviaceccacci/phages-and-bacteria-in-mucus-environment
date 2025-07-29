function [ numNod ] = giveNumNodesElementFromOrder(order,elementType)

    switch elementType

        case 'tri'   
            [ numNod ] = giveNumNodesElementFromOrder_tri(order) ; 

        case 'tet'
            [ numNod ] = giveNumNodesElementFromOrder_tet(order) ; 

        case 'quad'
            [ numNod ] = giveNumNodesElementFromOrder_quad(order) ; 
            
        otherwise
            [ numNod ] = (order + 1)^3;
            
    end
    
end