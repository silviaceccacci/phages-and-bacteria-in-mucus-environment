function [X2D,T2D] = generateRectangleMesh(ne1dx,ne1dy)

    x1d = linspace(-1,1,ne1dx+1);
    y1d = linspace(-1,1,ne1dy+1);
    
    np1dx = length(x1d);
    t1dx = [ 1:(np1dx-1) 
            2:np1dx ];
    np1dy = length(y1d);
    t1dy = [ 1:(np1dy-1) 
            2:np1dy ];
            
    X2D = zeros(2,np1dx*np1dy);
    T2D = zeros(4,ne1dx*ne1dy);
    for iy=1:np1dy
        if(iy<=ne1dy)
            list_e = (ne1dx*(iy-1)) + (1:ne1dx);
            T2D(1:2,list_e) = ((iy-1)*np1dx) + t1dx;
            T2D(3  ,list_e) = ((iy-0)*np1dx) + t1dx(2:2:end,:);
            T2D(4  ,list_e) = ((iy-0)*np1dx) + t1dx(1:2:end,:);
        end
        list_p = (np1dx*(iy-1)) + (1:np1dx);
        X2D(1,list_p) = x1d;
        X2D(2,list_p) = y1d(iy);
    end

    
    X2D = X2D';
    T2D = [T2D(1:3,:)  T2D([1,3,4],:)]';
 
end
