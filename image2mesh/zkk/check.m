function [] = check()
    close all;
    clear all;
    
%     addpath('./paraview')
    addpath('./input')
    addpath('./src/BAMG')
    addpath('./src/visualization')
    addpath('./src/grid')
    
    global machine;
    machine = 'abel';
    
    plotSteps = false;
    
    lambda = 5;
    

    fileName = 'airbos_f7_p5';
    format = 'jpg';
    
    imageName = [fileName '.' format];
    
    RGB = imread(imageName);
    
    figure(3)
    imshow(RGB)
    axis off
    axis equal
    colormap gray
    
    
    
	x0 =[ 0 0];
    
%     RGB = double(imread(imageName)); 
%     R = RGB(:,:,1);

    load('schlieren.mat', 'B');
    R = B;

%     R = double(rgb2gray(RGB));
    size(R)

    %error('para')
%      R = R(:, 1:200);
%      R = R(:, 2000:2500);
    
    %Z = 1.0*R' / 255.;
%     Z = 5*R' / 255;% / 100; 
    figure(10)
%     image(R)
%     colormap gray
   
    Z = R'/max(R(:));
    
    imshow(Z');
    
    Z = Z(:, end:-1:1);
           
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
            
    X = [x(:), y(:)];
    [T] = delaunay(x(:),y(:));
    
    %% Loop until no more elements are generated/deleted

    if(plotSteps)
        figure(1); clf; subplot(1,2,1); hold on
        plot(X(:,1),X(:,2),'r*')
        plotMesh(X,T,0)
    end

    M = [];

    scalar = lambda*Z(3:end-2,3:end-2);
    scalar = scalar(:);
    options.scalar = scalar;
    [X,T] = adapt(X,T,M,options);

    Xh = X;
    Th = T;
    %%
    
    figure(1)
    %plot(X(:,1),X(:,2),'*')
    %hold on
    plotMesh(Xh,Th,0)
    axis equal
    axis off
    %plot(Xh(:,1),Xh(:,2),'o')
%     title('Adapted mesh')

    figure(2)
    surface( x,y,Z(3:end-2,3:end-2),'EdgeColor','none')
    axis off
    axis equal
    colormap gray
%     title('Initial metric mesh')
    
    %X3D = [Xh, f(Xh(:,1),Xh(:,2))];
%     surf.z = Z(3:end-2,3:end-2);
%     surf.hx = hx;
%     surf.hy = hy;
%     surf.x0 = x0;
%     surf.width = width;
%     surf.height = height;
%     surf.structured = true;
%     X3D = [Xh,  160*project(surf, Xh, 'nodes')];
%     figure(3)
%     plotMesh(X3D,Th,0)
%     view(2)
%     title('Adapted "surface" mesh')
%     

%     mesh.X = X3D;
%     mesh.T = Th;
%     mesh.name = fileName;
%     exportTriMeshToParaview(mesh)    
%     
%     save(['./output/' fileName '.mat'],'mesh')
end
%%