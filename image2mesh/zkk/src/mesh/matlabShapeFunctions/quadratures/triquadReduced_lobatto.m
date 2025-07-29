%% 1r intent: no va be--> dona inf
% function [x,w] = triquadReduced_lobatto(N,vertices)
% 
%     [x,w] = lglnodesQuad(N);
%     index = find( (x(:,1)+x(:,2) )<=0);
%     x = x(index,:);
%     w = w(index);
% 
% %     figure(1)
% %     plot(x(:,1),x(:,2),'*')
% %     pause()
%     
%     tri1 = [  [-1;-1] [1;-1]  [-1;1]  ];
%     if(nargin<2)
%         tri2 = tri1;
%     else
%         tri2 = vertices';
%     end
%     [x2] = sendTri1PointsToTri2(x',tri1,tri2);
%     x = x2';
%       
% %     figure(2)
% %     plot(x(:,1),x(:,2),'*')
% %     pause()
%     
% end

%% 2n intent
function [x,w] = triquadReduced_lobatto(N,vertices)

%% Gauss-Lobatto on the square
    [x,w] = lglnodesQuad(N);   
    numPoints = length(w);
%     figure(1)
%     plot(x(:,1),x(:,2),'*')
%     pause()
%% Division of the points on two triangles
    index1 = find( (x(:,1)+x(:,2) )<=0);
    index2 = setdiff(1:numPoints,index1);
    tri1 = [  [-1;-1] [1;-1]  [-1;1]  ];
    x1 = x(index1,:);
    tri1bis = [ [-1;1] [1;-1]  [1;1] ];
    x2 = x(index2,:);
%% Mapping of the two triangles onto tri1end and tri1endbis    
    tri1end = [  tri1(:,2) [0;0] tri1(:,1)] ; %tri1end = [ tri1(:,1) tri1(:,2) [0;0] ] ;
    tri1endbis = [ tri1(:,1) [0;0] tri1(:,3) ] ;
    [x1bis] = sendTri1PointsToTri2(x1',tri1,tri1end);
%     figure(2)
%     plot(x1bis(1,:),x1bis(2,:),'*')
%     pause()  
    [x2bis] = sendTri1PointsToTri2(x2',tri1bis,tri1endbis);
%     figure(3)
%     plot(x2bis(1,:),x2bis(2,:),'*')
%     pause()
%% Composition of the two subtriangles: [-1,1]^2    
    x(index1,:) = x1bis';
    x(index2,:) = x2bis';
%% Mapping of the points on the reference triangle to the given one
    if(nargin<2)
        tri2 = tri1;
    else
        tri2 = vertices';
    end
    [x] = sendTri1PointsToTri2(x',tri1,tri2);
    x = x';
      
%     figure(4)
%     plot(x(:,1),x(:,2),'*')
%     pause()
% x 
% w
    
end