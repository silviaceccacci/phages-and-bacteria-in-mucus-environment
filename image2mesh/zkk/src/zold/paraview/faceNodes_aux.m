function res = faceNodes_aux(element)
    Tref = 1:element.numNod;
    switch element.type
        case 'tri'
            res = [...
                getEdge_tri(Tref,element.order,1) 
                getEdge_tri(Tref,element.order,2) 
                getEdge_tri(Tref,element.order,3) 
                ];
        case 'quad'
            res = [...
                getEdge_quad(Tref,element.order,1) 
                getEdge_quad(Tref,element.order,2) 
                getEdge_quad(Tref,element.order,3) 
                getEdge_quad(Tref,element.order,4) 
                ];
        otherwise
            error('Cannot plot this type of element yet')
    end
end