function p = orthopoly2D(x,n,element)
% p = orthopoly2D(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(xi,eta) in the reference triangle
%
% 
% xi = x(:,1); eta = x(:,2); 
% 
% % switch element.type
% %     case 'tri'
%         warning('OFF')
%         r = 2*(1+xi)./(1-eta)-1;
%         warning('ON')
%         % NOTE that to vectorize we have included this posteriori check, then the
%         % derived warning is not important
%         arrange=find(eta==1);  
%         r(arrange)=-1;
% 
%         s = eta;
%        
% %     case 'quad'
% %         r = xi;
% %         s = eta;
% % end

[r s] = xieta_to_rs(x(:,1),x(:,2));

p = orthopoly2D_rst([r,s],n,element);

end



% xi = x(:,1); eta = x(:,2); 
% 
% warning('OFF')
% r = 2*(1+xi)./(1-eta)-1;
% warning('ON')
% % NOTE that to vectorize we have included this posteriori check, then the
% % derived warning is not important
% arrange=find(eta==1);  
% r(arrange)=-1;
% 
% s = eta;
% 
% p = orthopoly2D_rst([r,s],n,element);
% 