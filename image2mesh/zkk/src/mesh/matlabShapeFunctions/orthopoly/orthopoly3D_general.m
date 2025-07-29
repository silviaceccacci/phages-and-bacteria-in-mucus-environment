function p = orthopoly3D_general(x,element)
    
    switch element.type
        case 'tet'
            p = orthopoly3D(x,element.order);
        case 'hex'
            p = orthopoly3D_hex(x,element.order);
    end
end


function p = orthopoly3D_hex(x,nDeg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x: points where we want to evaluade the shape functions
% element: element information of the approximation (order, distribution..)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [p_x]  = orthopoly1D(x(:,1), nDeg); 
    [p_y]  = orthopoly1D(x(:,2), nDeg);
    [p_z]  = orthopoly1D(x(:,3), nDeg);
   
    numNod = size(p_x,1)^3;
    p = zeros(numNod, size(p_x,2)); 
    count = 0;
    for ix = 1:size(p_x,1) %nDeg+1
        for iy = 1:size(p_y,1) %nDeg+1
            for iz = 1:size(p_z,1) %nDeg+1
                count = count+1;
                p(count,:)        =  p_x(ix,:).*   p_y(iy,:) .*  p_z(iz,:) ;
            end
        end
    end  
    
end

