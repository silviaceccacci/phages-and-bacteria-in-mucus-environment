function [] = generateMeshFromImage(fileName,format)

    plotSteps = false;
    
    global c;
    c = 36;%00000;
    
%     fileName = 'airbos_f7_p5';
%     format = 'jpg';
    
    imageName = [fileName '.' format];
    
    %load image
%     x0 =[ -2*pi -2*pi];
%     width = 4*pi;
%     height = 4*pi;
    
%     x0 =[ -1 -1];
%     width = 2;
%     height = 2;
    

%    m = 256;
%    n = m;

	x0 =[ 0 0];
    
    %RGB = double(imread('marmousi_hard.gif'));
    RGB = double(imread(imageName)); 
    R = RGB(:,:,1);
    size(R)
    %error('para')
    %R = R(1:1990,1:1990);
    %R = R(1200:1700,700:1200);
    
    %Z = 1.0*R' / 255.;
%     Z = 5*R' / 255;% / 100; 
    figure(10)
%     image(R)
%     colormap gray
   
    Z = R'/255;
    
    imshow(Z');
    
    Z = Z(:, end:-1:1);
        
    %Z = Z(200:500,100:400); 
    %Z = Z(1:500,100:400); 
    %Z = Z(1:end,100:400); 
    
    %Z = Z(500:1500,500:1500); 
    
    [m, n, channels] = size(Z);
    
    width  = m;
    height = n;  

    [x,y] = meshgrid(0:(m-1),0:(n-1));
    x = x';
    y = y';
    %[x,y] = meshgrid(0:(n-1),0:(m-1));
           
    hx = width/m;
    hy = height/n;    
    
    x = hx*x + x0(1);
    y = hy*y + x0(2);
    
    x = x(3:end-2,3:end-2);
    y = y(3:end-2,3:end-2);
    
    x0 = [ min(x(:)) min(y(:)) ];
    
    %[Z] = f(x,y);
    
    %[M0]  =  metric(Z,hx,hy);
    H0 = hessian(Z, hx , hy);
    M0 = metricInPointsFromHessianInPoints(H0);
    
    pixelPoints =[ (x(1:end-1) + x(2:end))/2.0  ;  (y(1:end-1) + y(2:end))/2.0 ]';
    
    field0.z = M0;
    field0.hx = hx;
    field0.hy = hy;
    field0.x0 = x0;
    field0.width = width;
    field0.height = height;
    field0.x = x;
    field0.y = y;
    field0.pixelPoints = pixelPoints;
    
    field0  
    
    hessianField0 = field0;
    hessianField0.z = H0;
              
    hessianField0
        
    X = [x(:), y(:)];
    [T] = delaunay(x(:),y(:));
    
%     nx = m-1;
%     ny = n-1;
%     numQuads = nx*ny;
%     l = zeros(m,n);
%     l(:) = 1:(m*n);
%     
%     l1 = l(1:end-1,1:end-1);
%     l2 = l(2:end,1:end-1);
%     l3 = l(2:end,2:end);
%     l4 = l(1:end-1,2:end);
%     
%     tquad = zeros(numQuads,4);    
%     tquad(:,:) = [l1(:) l2(:) l3(:) l4(:)];    
% %     ielem = 0;
% %     for j=1:ny
% %         for i=1:nx
% %             ielem = ielem +1;
% %             tquad(ielem,:) = [l(i,j) l(i+1,j) l(i+1,j+1) l(i,j+1)];
% %         end
% %     end
%     T = [tquad(:,[1 2 4]); tquad(:,[ 2 3 4])];
    %% Loop until no more elements are generated/deleted
    maxIt = 1;
    
    it=0;
    converged = false;
    numElemsOld = size(T,1);
    while(it<maxIt && converged==false)
    %while(it<maxIt)
    
        if(plotSteps)
            figure(1); clf; subplot(1,2,1); hold on
            plot(X(:,1),X(:,2),'r*')
            plotMesh(X,T,0)
        end
        
        % %M = Metric(@(x,y)f(x,y), X(:,1),X(:,2));
        
        M = [];
%         if(it==0)
%             disp('Computing metric')
%             M = project(field0, X);
%         else
%             %M = smartProject(field0,X,T);
%             H = smartProject(hessianField0,X,T);
%             disp('Building metric from smart projected hessian')
%             M = metricInPointsFromHessianInPoints(H);
%         end
        
        if(it==0)
        scalar = 160*Z(3:end-2,3:end-2);
        scalar = scalar(:);
        options.scalar = scalar;
        [X,T] = adapt(X,T,M,options);
        else
%         surf.z = Z(3:end-2,3:end-2);
%         surf.hx = hx;
%         surf.hy = hy;
%         surf.x0 = x0;
%         surf.width = width;
%         surf.height = height;
%         scalar = project(surf, Xh);
%         scalar = scalar(:);
%         options.scalar = scalar;
%         [X,T] = adapt(X,T,M,options);           
        end

        it = it+1;
        numElems = size(T,1);
        converged = abs(numElemsOld-numElems)/numElems < 0.001;%(numElemsOld + numElemsOld*0.05) > numElems;
        numElemsOld = numElems;
        
        if(plotSteps)
            figure(1); subplot(1,2,2); hold on
            plot(X(:,1),X(:,2),'r*')
            plotMesh(X,T,0)
            pause(2)
        end
        
        
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
    surface( x,y,Z(3:end-2,3:end-2),'EdgeColor','none')
    title('Initial metric mesh')
    
    %X3D = [Xh, f(Xh(:,1),Xh(:,2))];
    surf.z = Z(3:end-2,3:end-2);
    surf.hx = hx;
    surf.hy = hy;
    surf.x0 = x0;
    surf.width = width;
    surf.height = height;
    X3D = [Xh,  160*project(surf, Xh)];
    figure(3)
    plotMesh(X3D,Th,0)
    view(2)
    title('Adapted "surface" mesh')
    

    mesh.X = X3D;
    mesh.T = Th;
    mesh.name = fileName;
    exportTriMeshToParaview(mesh)    
    
    save(['./output/' fileName '.mat'],'mesh')
end
%%
function [Z] = f(X,Y)
    %Z = X.^2 + Y.^2;
    
    %Z = cos(4.*sqrt((X+0.1).^2 + (Y+0.1).^2));
    %Z = cos(1.*sqrt((abs(X)+0.1).^2 + (abs(Y)+0.1).^2));
    
    %Z = 100*X.^2 + Y.^2;
    
    %a = sqrt(2)/2.0;
    %Z = 10*(a*X - a*Y).^2 + (a*X + a*Y).^2;
    
    %a = 0.5; 
    %b = 0.1;
    %Z = sqrt((X.^2)/(a.^2) + (Y.^2)/(b.^2));
    
    
    %Z = (10*X.^3 + Y.^3) + atan2(0.01,(sin(5*Y) - 2*X));
    
    
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

    H = hessian(f, hx , hy);

    [M] = metricInPointsFromHessianInPoints(H);
    
%     M = zeros(size(H));
%     for j=1:size(H,2)
%         for i=1:size(H,1)
%            A =  [H(i,j,1) H(i,j,2); H(i,j,2) H(i,j,3)];
% 
%            M(i,j,:) = metricFromHessian(A);
%         end
%     end  
       
    
%     global c;
%     H = hessian(f, hx , hy);
% 
%     for j=1:size(H,2)
%         for i=1:size(H,1)
%            A =  [H(i,j,1) H(i,j,2); H(i,j,2) H(i,j,3)];
%            [V,D] = eig(A);
%            A = V*abs(D)*V';
%            H(i,j,1) = A(1,1);
%            H(i,j,2) = A(2,1);
%            H(i,j,3) = A(2,2);
%         end
%     end  
%        
%     M = c .*(H); 
end

function [M] = metricFromHessian(H)
    global c;
    
    M = zeros(3,1);
    
    [V,D] = eig(H);
    
    H = V*abs(D)*V';
    
    M(1) = H(1,1);
    M(2) = H(2,1);
    M(3) = H(2,2);
    
    M = c*M;
end

function [M] = metricInPointsFromHessianInPoints(H)


    M = zeros(size(H));

    if(ndims(H)==3)
        for j=1:size(H,2)
            for i=1:size(H,1)
               A =  [H(i,j,1) H(i,j,2); H(i,j,2) H(i,j,3)];
               M(i,j,:) = metricFromHessian(A);
            end
        end  
    elseif(ndims(H)==2)
        for i=1:size(H,1)
           A =  [H(i,1) H(i,2); H(i,2) H(i,3)];
           M(i,:) = metricFromHessian(A);
        end          
    else
        error('not recognised hessian dimensions');
    end

end














