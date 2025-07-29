function [p,dp_dxi,dp_deta,dp_dzeta] = orthopoly3D_deriv_xietazeta(x,n)
%
% [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

% AIXO MHO HE PATILLAT JO

% xi = x(:,1); eta = x(:,2); zeta = x(:,3);
% 
% r = - ones(size(xi));
% s =   ones(size(xi));
% 
% ind_aux    = find( (eta+zeta)~=0  )  ;
% ind3       = ind_aux( find( zeta(ind_aux) ~= 1 ) ) ;
% 
% r(ind3) = -2*(1+ xi(ind3))./(eta(ind3)+zeta(ind3))-1;
% s(ind3) =  2*(1+eta(ind3))./(1-zeta(ind3))-1;
% 
% t = zeta;

[r s t] = xietazeta_to_rst(x(:,1), x(:,2),x(:,3)); % canviar nom funcio

[p,dp_dxi,dp_deta,dp_dzeta] = orthopoly3D_deriv_rst([r s t],n);

end


