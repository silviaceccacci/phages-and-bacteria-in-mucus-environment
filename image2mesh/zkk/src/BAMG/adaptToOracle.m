function [] = adaptToImage()
    close all;
    clear all;
    
    global c;
    c = 25;%00000;
    
    %load image
    m = 100;
    n = m;
    %[x,y] = meshgrid(-2:(m+2),-2:(n+2));
    [x,y] = meshgrid(-(m+2):(m+2),-(n+2):(n+2));
    %[x,y] = meshgrid(-(2):(m+2),-(2):(n+2));
    
    %extend image!!!
    
    hx = 2*pi/m;
    hy = 2*pi/n;
    
    %hx = 1.0/m;
    %hy = 1.0/n;    
    
    x = hx*x;
    y = hy*y;
        
    [Z] = f(x,y);
    [M]  =  metric(Z,hx,hy);
    
    x = x(3:end-2,3:end-2);
    y = y(3:end-2,3:end-2);
    
    X = [x(:), y(:)];
    [T] = delaunay(x(:),y(:));
        
    %% One iteration
    % M = reshape(M,size(M,1)*size(M,2),size(M,3));
    %[Xh,Th] = adapt(X,T,M);
   
% %     Xini = X;
% %     Tini = T;
% %     mesh.X = Xini;
% %     mesh.T = Tini;
% %     options.meshToAdapt = mesh;
% %     [Xh,Th] = adapt(X,T,M,options);

    %% Loop until no more elements are generated/deleted
    maxIt = 100;
    
    it=0;
    converged = false;
    numElemsOld = size(T,1);
    while(it>maxIt || converged==false)
    %while(it<maxIt)
        
        figure(1); clf; subplot(1,2,1); hold on
        plot(X(:,1),X(:,2),'r*')
        plotMesh(X,T,0)
    
        M = Metric(@(x,y)f(x,y), X(:,1),X(:,2));
        [X,T] = adapt(X,T,M);

        it = it+1;
        numElems = size(T,1);
        converged = abs(numElemsOld-numElems)/numElems < 0.001;%(numElemsOld + numElemsOld*0.05) > numElems;
        numElemsOld = numElems;
        
        figure(1); subplot(1,2,2); hold on
        plot(X(:,1),X(:,2),'r*')
        plotMesh(X,T,0)
        pause(2)
        
        
        [it numElems]
    end
    Xh = X;
    Th = T;
    %%
    
    figure(1)
    %plot(X(:,1),X(:,2),'*')
    %hold on
    plotMesh(Xh,Th,0)
    %plot(Xh(:,1),Xh(:,2),'o')
    title('Adapted mesh')

    figure(2)
    surface( x,y,Z(3:end-2,3:end-2))%,'EdgeColor','none')
    title('Initial metric mesh')
    
    X3D = [Xh, f(Xh(:,1),Xh(:,2))];
    figure(3)
    plotMesh(X3D,Th,0)
    view(2)
    title('Adapted "surface" mesh')
    
end
%%
function [Z] = f(X,Y)
    %Z = X.^2 + Y.^2;
    
    %Z = cos(4.*sqrt((X+0.1).^2 + (Y+0.1).^2));
    Z = cos(1.*sqrt((abs(X)+0.1).^2 + (abs(Y)+0.1).^2));
    
    %Z = 100*X.^2 + Y.^2;
    
    %a = sqrt(2)/2.0;
    %Z = 10*(a*X - a*Y).^2 + (a*X + a*Y).^2;
    
    %a = 0.5; 
    %b = 0.1;
    %Z = sqrt((X.^2)/(a.^2) + (Y.^2)/(b.^2));
    
    %Z = (10*X.^3 + Y.^3) + atan2(0.01,(sin(5*Y) - 2*X));
end

%% continuous
function [fx] = Dx(f, x, y)
    h = 10.^(-16/3.0);
    fx = (f(x + h, y) - f(x - h, y)) ./ (2*h);
end 

function [fy] = Dy(f, x, y)
    h = 10^(-16/3.0);
    fy = (f(x, y + h) - f(x, y - h)) ./ (2*h);
end

function [hessf] = Hessian(f, x, y)
    fxx = Dx(@(x,y)Dx(f, x, y), x, y);
    fxy = Dy(@(x,y)Dx(f, x, y), x, y);
    fyy = Dy(@(x,y)Dy(f, x, y), x, y);
    
%     h = 0.001;%(10^(-16/9.0));
%     fxx = (f(x+2*h,y) - 2*f(x,y) + f(x-2*h,y) ) ./ (4*h^2);
%     fyy = (f(x,y+2*h) - 2*f(x,y) + f(x,y-2*h) ) ./ (4*h^2);
%     fxy = (f(x+h,y+h) - f(x-h,y+h) - f(x+h,y-h) + f(x-h,y-h)) ./ (4*h^2);
    
    hessf = cat(2,fxx,fxy,fyy);
end

function [M]  =  Metric(f, x, y)
    global c;
    H = Hessian(f, x, y);
    
    for i=1:size(x,1)
       A =  [H(i,1) H(i,2); H(i,2) H(i,3)];
       [V,D] = eig(A);
       A = V*abs(D)*V';
       H(i,1) = A(1,1);
       H(i,2) = A(2,1);
       H(i,3) = A(2,2);
    end  
        
    M = c .* (H);
end

%% discrete
function [fx] = dx(f, h)       
    fnext = f(3:end,:);
    fprev = f(1:end-2,:);
    
    fx = (fnext - fprev) / (2*h);
end

function [fy] = dy(f, h)       
    fnext = f(:,3:end);
    fprev = f(:,1:end-2);
    
    fy = (fnext - fprev) / (2*h);
end

function [hessf] = hessian(f, hx , hy)
    fxx = dx(dx(f,hx),hx);
    fyy = dy(dy(f,hy),hy);
    fxy = dy(dx(f,hx),hy);
    
    fxx = fxx(:, 3:end-2);
    fyy = fyy(3:end-2, :);
    fxy = fxy(2:end-1, 2:end-1);
    
    hessf = cat(3,fxx,fxy,fyy);
end

function [M]  =  metric(f, hx , hy)
    global c;
    H = hessian(f, hx , hy);

    for j=1:size(H,2)
        for i=1:size(H,1)
           A =  [H(i,j,1) H(i,j,2); H(i,j,2) H(i,j,3)];
           [V,D] = eig(A);
           A = V*abs(D)*V';
           H(i,j,1) = A(1,1);
           H(i,j,2) = A(2,1);
           H(i,j,3) = A(2,2);
        end
    end  
       
    M = c .*(H); 
end


%% grid finding
function [location] = locate(x,field)
    
    hx = field.hx;
    hy = field.hy;
    z  = field.z;
    [m n] = size(z);
    
    i = floor(x(:,1)/hx);
    j = floor(x(:,2)/hy);
    
    location.l = i+(j-1)*m;
    
    %chi = [x(:,1)/hx  - i ,  x(:,2)/hy  - j] ;
    
end


function [z_on_x] = project(field, x)    
    [location] = locate(x, field);
    
    z = field.z;
    [m, n, numComponents] = size(z);
    
    numPoints = size(x, 1);
    
    z_on_x = zeros(numPoints, numComponents);
    
    for c = 1:numComponents
        zc = z(:,:,c);        
        z_on_x(:, c) = zc(location.l);
    end    
end











