function [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n,element)
%
% [p,dp_dxi,dp_deta] = orthopoly2D_deriv_rst(x,n)
% Computes the ortogonal base of 2D polynomials of degree less 
% or equal to n at the point x=(r,s) in [-1,1]^2
%

% N = (n+1)*(n+2)/2 ;%number of nodes/polynomials
N = element.numNod;
p = zeros(N,size(x,1));
dp_dxi  = zeros(N,size(x,1));
dp_deta = zeros(N,size(x,1));

r = x(:,1); s = x(:,2); 

xi = (1+r).*(1-s)/2-1;
eta = s;

% facTol  = 1e-14;
% arrange=find(eta==1);
% s(arrange)   = s(arrange)-facTol;
% eta(arrange) = eta(arrange)-facTol;

if(find(eta==1))
   s = s*0.9999999999; 
   eta = eta*0.9999999999;
end

% warning('OFF')
dr_dxi  = 2./(1-eta);
dr_deta = 2*(1+xi)./(1-eta).^2;
% warning('ON')
% disp('canvi a orthopoly2D_deriv_rst linies 23-25')
% arrange=find(eta==1);
% % dr_dxi(arrange)=0;%-1;%0;%-1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PER QUIN VALOR SHA DE CANVIAR?
% % dr_deta(arrange)=0;%1;%0;%1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PER QUIN VALOR SHA DE CANVIAR?
% facTol  = 1e-12;
% dr_dxi(arrange)  = 2./(1-eta(arrange)+facTol);
% dr_deta(arrange) = 2*(1+xi(arrange))./(1-eta(arrange)+facTol).^2;


%Ordering: 1st incresing the degree and 2nd lexicogafic order
ncount = 0; %counter for the polynomials order
%Loop on degree
for nDeg = 0:n
  %Loop increasing i
  for i = 0:nDeg
     if i==0
         p_i = ones(size(x,1),1); %1;
         q_i = ones(size(x,1),1); % 1;
         dp_i = zeros(size(x,1),1); %0;
         dq_i = zeros(size(x,1),1); %0;
     else
         p_i = jacobiP_vect(r,0,0,i);
         dp_i = jacobiP_vect(r,1,1,i-1)*(i+1)/2;    
%          dq_i = q_i.*(-i)/2;
         q_i = q_i.*(1-s)./2;
         dq_i = q_i.*(-i)./(1-s); %introduint qi en dqi es pot calcular millor
% % %          disp('canvi a orthopoly2D_deriv_rst linies 44-45')
%          arrange = find(s==1);
% %          dq_i(arrange)=0;%-1; %%%%%%%%%%%%%%%%%%%%%%%%%%%%% PER QUIN VALOR SHA DE CANVIAR?
%         dq_i(arrange) = q_i(arrange).*(-i)./(1-s(arrange)-facTol);
     end
     %Value for j
     j = nDeg-i;
     if j==0
        p_j = ones(size(x,1),1); %1;
        dp_j = zeros(size(x,1),1); %0; 
     else
        p_j = jacobiP_vect(s,2*i+1,0,j);
        dp_j = jacobiP_vect(s,2*i+2,1,j-1)*(j+2*i+2)/2;  
     end
     ncount= ncount+1;
     factor = sqrt( (2*i+1)*(i+j+1)/2 );
     %Normalized polinomial
     p(ncount,:)    = ( p_i.*q_i.*p_j )*factor;
     %Derivatives with respect to (r,s)
     dp_dr = ( (dp_i).*q_i.*p_j )*factor; 
     dp_ds = ( p_i.*(dq_i.*p_j+q_i.*dp_j) )*factor;
     %Derivatives with respect to (xi,eta)
     dp_dxi(ncount,:)  = dp_dr.*dr_dxi;
     dp_deta(ncount,:) = dp_dr.*dr_deta + dp_ds;
  end
end

