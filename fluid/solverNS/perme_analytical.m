function [kval]=perme_analytical(X)

kx=por1D(X(1));
ky=por1D(X(2));
kval = kx*ky;

% kval = 1+kval;
% kval = 1/kval;
% kval = (kval -0.5)*2.0;
% kval = kval^0.3;
% kmin = 0.01;
% kmax = 1+kmin;
% kval = (kval+kmin)/kmax;

% kval=kval*1000;
% kval=kval+1e-2;
% kval = kval^0.5;
% kval = kval+1e-6;
% kval = 1/kval;
kval = kval^20;
kval = kval*1000;
kmin = 1e-5;
kval = kval+kmin;
kval = 1/kval;

% kval=kval*100;
% kval=1+kval;

% xplot = linspace(-1,1,100);
% yplot = por1D(xplot);
% figure;
% plot(xplot,yplot)
% pause()


% A = [sin(pi*X(1))    0
%          0            sin(pi*X(2))]; % create funciton K(X) interpolating in mesh
% 
% A = (A+1)/2;
% kmin = 0.5;
% kmax = 1+kmin;
% A = (A+kmin)/kmax;
end

function kval=por1D(x)
    %kval = sin(pi*X(1)*X(2));
    kval = cos(pi*x/2.0);%*sin(pi*X(2)/2.0);
    %kval = (kval+1)/2;
    %kval = abs(kval);
    
%     kval = 2/(1+kval);
% 
%     kmin = 0.5;
%     kmax = 1+kmin;
%     kval = (kval+kmin)/kmax;

end