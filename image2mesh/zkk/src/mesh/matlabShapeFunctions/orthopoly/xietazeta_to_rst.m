function [a,b,c] = xietazeta_to_rst(r,s,t) 
    
    Np = length(r);
    a = zeros(Np,1); b = zeros(Np,1);
    for n=1:Np  
      if(s(n)+t(n) ~= 0)
        a(n) = 2*(1+r(n))/(-s(n)-t(n))-1;
      else
        a(n) = -1;
      end
      if(t(n) ~= 1)
        b(n) = 2*(1+s(n))/(1-t(n))-1;
      else
        b(n) = -1;
      end
    end
    c = t;
    
%     problemScale = 0.9999;%0.999999999999999;
    problemScale = 0.9999999999999;%0.999999999999999;
% %     problemScale = 0.95;%0.999999999999999;

    a = a*problemScale ;
    b = b*problemScale ;
    c = c*problemScale ;
    
    
        
%     problemScale = 0.999999999999999;
%     Np = length(r);
%     a = zeros(Np,1); b = zeros(Np,1);
%     for n=1:Np  
%       if(s(n)+t(n) ~= 0)
%         a(n) = 2*(1+r(n))/(-s(n)-t(n))-1;
%       else
%         a(n) = -1;
%       end
%       if(t(n) ~= 1)
%         b(n) = 2*(1+s(n))/(1-t(n))-1;
%       else
%         b(n) = -1;
%       end
%     end
%     c = t;
%     ind=find(abs(a)==1);
%     a(ind) = a(ind)*problemScale;
%     ind=find(abs(b)==1);
%     b(ind) = b(ind)*problemScale;
%     ind=find(abs(c)==1);
%     c(ind) = c(ind)*problemScale;
    
end