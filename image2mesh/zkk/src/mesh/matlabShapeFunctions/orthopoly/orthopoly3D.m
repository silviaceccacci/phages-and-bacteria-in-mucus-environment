function p = orthopoly3D(x,n)
% p = orthopoly3D(x,n)
% Computes the ortogonal base of 3D polynomials of degree less 
% or equal to n at the point x=(xi,eta,zeta) in the reference tetrahedra
%

% xi = x(:,1); eta = x(:,2); zeta = x(:,3);
% 
% r = - ones(size(xi));
% s =   ones(size(xi));
% 
% ind_aux    = find( (eta+zeta)~=0  )  ;
% ind3       = ind_aux( find( zeta(ind_aux) ~= 1 ) ) ;
% 
% r(ind3) = -2*(1+ xi(ind3))./(eta(ind3)+zeta(ind3)) - 1;
% s(ind3) =  2*(1+eta(ind3))./(1-zeta(ind3)) - 1;
% 
% t = zeta;

[r s t] = xietazeta_to_rst(x(:,1), x(:,2),x(:,3)); % canviar nom funcio

p = orthopoly3D_rst([r,s,t],n);

if(find(isnan(p)))
    error('mierda')
end
end

%% checkejar a veure si hi he ficat la pota
% xi = x(1); eta = x(2); zeta = x(3);
% 
% if (eta+zeta)==0 
%     r = -1; s=1;
% elseif zeta==1
%     r = -1; s=1;  %or s=-1 (check that nothing changes)
% else
%     r = -2*(1+xi)/(eta+zeta)-1;
%     s = 2*(1+eta)/(1-zeta)-1;
% end
% t = zeta;
% 
% p = orthopoly3D_rst([r,s,t],n);


%% Code with comments
% xi = x(:,1); eta = x(:,2); zeta = x(:,3);
% 
% r = - ones(size(xi));
% s =   ones(size(xi));
% 
% % ind1       = find( (eta+zeta)==0 ) ) ;
% ind_aux    = find( (eta+zeta)~=0  )  ;
% % ind2       = ind_aux( find( zeta(ind_aux) == 1 ) ) ;
% ind3       = ind_aux( find( zeta(ind_aux) ~= 1 ) ) ;
% 
% % r(ind2) = -1;
% % s(ind2) = -1;
% 
% r(ind3) = -2*(1+ xi(ind3))./(eta(ind3)+zeta(ind3))-1;
% s(ind3) =  2*(1+eta(ind3))./(1-zeta(ind3))-1;
% 
% t = zeta;
% 
% p = orthopoly3D_rst([r,s,t],n);
% 
% % p = orthopoly3D_rst(x,n);
