function [e] = getEdge(T, arg2, localEdge)
% arg2 -> element    /    order


    if(isstruct(arg2))
        element = arg2;
        order = element.order;
    else
        order = arg2;
        numNodesElement = size(T,2);
        if(     numNodesElement == giveNumNodesElementFromOrder_tri(order)  )
            element.type = 'tri' ;   
        elseif( numNodesElement == giveNumNodesElementFromOrder_tet(order)  ) 
            element.type = 'tet';
        elseif( numNodesElement == giveNumNodesElementFromOrder_quad(order)  ) 
            element.type = 'quad';
        elseif( numNodesElement == giveNumNodesElementFromOrder_hex(order)  ) 
            element.type = 'hex';
        end
    end


    switch  element.type

        case 'tri'
            [e] = getEdge_tri(T, order, localEdge);

        case 'tet'
            [e] = getEdge_tet(T, order, localEdge);

        case 'quad'
            [e] = getEdge_quad(T, order, localEdge);
            
        case 'hex'
            [e] = getEdge_hex(T, order, localEdge);
            
        otherwise
            error('Something is wrong: sure element selected or connectivity correct?')
    end
    
end