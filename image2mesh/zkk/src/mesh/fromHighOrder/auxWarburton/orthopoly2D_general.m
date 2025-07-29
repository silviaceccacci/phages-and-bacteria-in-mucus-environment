function p = orthopoly2D_general(x,n,element)

    switch element.type
        case 'tri'
            p = orthopoly2D(x,n,element);
        case 'quad'
            p = orthopoly2D_quad(x,n,element);
    end

end