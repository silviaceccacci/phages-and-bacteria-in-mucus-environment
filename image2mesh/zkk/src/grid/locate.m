function [location] = locate(x,field)
    
    hx = field.hx;
    hy = field.hy;
    x0 = field.x0;
    if(isfield(field,'z'))
        [m n k] = size(field.z);
    else
        m = field.m;
        n = field.n;
    end
    
    i = floor((x(:,1)-x0(:,1))/hx)+1;
    j = floor((x(:,2)-x0(:,2))/hy)+1;
     
%     i = mod(i-1,m)+1;
%     j = mod(j-1,n)+1;

    %if(max(i<1) || max(j<1) || max(i>m) || max(j>n) )
    if(~isempty(find(i<1)) || ~isempty(find(j<1)) || ~isempty(find(i>m)) || ~isempty(find(j>n)) )
        error('Out of the field data');
    end
    
    location.l = i + (j-1)*m;
    location.i = i ;
    location.j = j;

    %location.l = j + (i-1)*m;
    
    %location.l = j + (i-1)*n;
    %location.l = j + (i-1)*n;
    %location.l = i + (j-1)*n;
    
    %chi = [x(:,1)/hx  - i ,  x(:,2)/hy  - j] ;
    
end