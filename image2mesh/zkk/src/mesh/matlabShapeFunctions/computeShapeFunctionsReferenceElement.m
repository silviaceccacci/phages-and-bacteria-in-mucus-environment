function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement(nDeg,coord,option,element)
%
% [shapeFunctions,gaussWeights]=computeShapeFunctionsReferenceElement(nDeg,
% coord)
%
% Function to compute the shape functions (& derivatives) at Gauss points
%
% Input:
% nDeg:  degree of interpolation
% coord: nodal coordinates at the reference element
%
% Output:
% shapeFunctions: shape functions evaluated at the gauss points
%                 size is nOfNodes X nOfGauss X (nsd + 1)
%                 nsd+1 because shape function (1)
%                 and derivatives (nsd) are stored
% gaussWeights:   weights
%

nsd = size(coord,2);
if nsd==1
    nOfGaussPoints = option;
    [shapeFunctions,gaussWeights,gaussPoints]=...
        computeShapeFunctionsReferenceElement1D(nDeg,coord,nOfGaussPoints);
elseif nsd==2
    coordSubElement = option;
    [shapeFunctions,gaussWeights,gaussPoints]=...
        computeShapeFunctionsReferenceElement2D(nDeg,coord,coordSubElement,element);
elseif nsd==3
    [shapeFunctions,gaussWeights,gaussPoints]=...
        computeShapeFunctionsReferenceElement3D(nDeg,coord,option,element);
else
    error('wrong nsd in computeShapeFunctionsReferenceElement')
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement1D(nDeg,coord,nOfGaussPoints)

%number of nodes/polynomials
nOfNodes = nDeg+1;
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement1D')
end

% [z,w] = gaussLegendre(nOfGaussPoints,-1,1);
[z,w] = lglnodes(nOfGaussPoints-1);
nOfGauss = length(w);

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,2);
gaussWeights = w';
% gaussPoints = zeros(nOfGauss,1);

%Integration over [-1,1]
% x = z' ;
x =z;
[p,p_xi] = orthopoly1D_deriv(x,nDeg);
% N = U\(L\(P*[p,p_xi]));
shapeFunctions(:,:,1) = U\(L\(P*p));
shapeFunctions(:,:,2) = U\(L\(P*p_xi));
% only for PFEM
gaussPoints = x';

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement2D(nDeg,coord,coordElem,element)
%% TRIANGLES and QUADRILATERALS
    switch element.type
    case 'tri'
        if( length(coordElem)<=1 || isempty( find( coordElem~= [-1 -1; 1 -1; -1 1])  ) )
            vert = [-1 -1; 1 -1; -1 1];
        end
    case 'quad'
        if( length(coordElem)<=1 || isempty( find( coordElem~= [-1 -1; 1 -1; 1 1 ;-1 1])  ) )
            vert = [-1 -1; 1 -1; 1 1 ;-1 1];
        end
    end
    if(isfield(element,'orderQuadrature'))
        N = element.orderQuadrature;
    else
        N = element.order+1;%= 2*nDeg+1; 
    end

    % quadrature
    element.vert = vert;
    [x,w] = getQuadrature(element,N);
    
% figure(13)
%     plot(x(:,1),x(:,2),'*')
%     hold on
%     plot(coord(:,1),coord(:,2),'o')
%     pause()
    % get shape functions
    [shapeFunctions]=getShapeFunctions(element,coord,x);
    gaussWeights = w;
    gaussPoints = x;

    if(find(isnan(shapeFunctions)))
        error('is nan in shapeFun 2D')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement3D(nDeg,coord,options,element)
%% Righ now just tetrahedra, but for hexes everything is already done
if(isfield(options,'vert')==false || max( isempty(options.vert) ) ||  (options.vert(1,1) == 0) || ...
        isempty( find( options.vert~= [ 1 -1 -1;-1 -1 -1;-1  1 -1;-1 -1  1 ])  ) )
    vert = [ 1 -1 -1
            -1 -1 -1
            -1  1 -1
            -1 -1  1 ];
        
        if(isfield(element,'orderQuadrature'))
            N = element.orderQuadrature;
        else
%             N = 2*nDeg+1; %     N = nDeg;
            N = element.order+1;%ceil((nDeg+1)/2); %nDeg+1;%nDeg+1;
        end
else
    vert = options.vert;
    N = 2*ceil(sqrt(nDeg)) +1;
end

% quadrature
element.vert = vert;
[x,w] = getQuadrature(element,N);

% get shape functions
[shapeFunctions]=getShapeFunctions(element,coord,x);
gaussWeights = w;
gaussPoints = x;

if(find(isnan(shapeFunctions)))
    error('is nan in shapeFun 3D')
end

% figure(100)
% plot3(gaussPoints(:,1),gaussPoints(:,2),gaussPoints(:,3),'*r')
% pause()

end





%% Possibilitats triangles:
%% una mica mes antiga q lactual
% 
% %number of nodes/polynomials
% nOfNodes = (nDeg+1)*(nDeg+2)/2;
% if nOfNodes~=size(coord,1)
%     error('Error computeShapeFunctionsReferenceElement2D')
% end
% 
% if( length(coordTri)<=1 || isempty( find( coordTri~= [-1 -1; 1 -1; -1 1])  ) )  %  w = 2*w; z = 2*z -1; %mapping onto the normal reference triangle
%     N=nDeg^2;                    %4*nDeg*nDeg;%nDeg; %--> sembla bastar amb 2*nDeg
%     [X,Y,Wx,Wy]=triquad(N,[-1 -1; 1 -1; -1 1]);
%     z=zeros(N*N,2);
%     w=zeros(1,N*N);
%     count=0;
%     for i=1:N
%         for j=1:N
%             count=count+1;
%             z(count,:)=[X(i,j) Y(i,j)];
%             w(count)=Wx(i)*Wy(j);
%         end
%     end  
% else
%     % THIS CODE LETS US SELECT THE DESIRED NUMBER OF INTEGRATION POINTS
%     N=nDeg;                %2*nDeg;%fix(sqrt(nDeg));
%     [X,Y,Wx,Wy]=triquad(N,[coordTri(1,:); coordTri(2,:); coordTri(3,:) ]);
%     z=zeros(N*N,2);
%     w=zeros(1,N*N);
%     count=0;
%     for i=1:N
%         for j=1:N
%             count=count+1;
%             z(count,:)=[X(i,j) Y(i,j)];
%             w(count)=Wx(i)*Wy(j);
%         end
%     end  
% end
% 
% nIP = length(w); %number of integration points
% nOfGauss = nIP; % for the cubature
% %Vandermonde matrix
% V = Vandermonde_LP(nDeg,coord);
% [L,U,P] = lu(V');
% 
% shapeFunctions = zeros(nOfNodes,nOfGauss,3);
% x = z; % (r,s) coordinates
% [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg);
% shapeFunctions(:,:,1) = U\(L\(P*p));
% shapeFunctions(:,:,2) = U\(L\(P*p_xi));
% shapeFunctions(:,:,3) = U\(L\(P*p_eta));
% gaussWeights = w';
% gaussPoints = x;
% % plot(x(:,1),x(:,2),'*')
% % pause()
%% Funcio tabulada tris (david)

% %number of nodes/polynomials
%     nOfNodes = (nDeg+1)*(nDeg+2)/2;
%     if nOfNodes~=size(coord,1)
%         error('Error computeShapeFunctionsReferenceElement2D')
%     end
%    switch nDeg
%         case 1
%             OrderCubature = 5; % 7 pdG
%         case 2
%             OrderCubature = 10;
%         case 3
%             OrderCubature = 10;
%         case 4
%             OrderCubature = 15;
%         case 5
%             OrderCubature = 15;
%         case 6
%             OrderCubature = 15;
%         case 7
%             OrderCubature = 15;
%         case {8,9,10,11}
%             OrderCubature = 25;
%     end
% 
%     if( length(coordTri)<=1 )
% %         THIS CODE DOES NOT LET US CHOOSING A LOT OF INTEGRATION POINTS
%         [z,w] = GaussLegendreCubature2D(OrderCubature);
%         w = 2*w; z = 2*z -1; %mapping onto the normal reference triangle
%         N=nDeg^2;                    %4*nDeg*nDeg;%nDeg; %--> sembla bastar amb 2*nDeg
%         [X,Y,Wx,Wy]=triquad(N,[-1 -1; 1 -1; -1 1]);
%         z=zeros(N*N,2);
%         w=zeros(1,N*N);
%         count=0;
%         for i=1:N
%             for j=1:N
%                 count=count+1;
%                 z(count,:)=[X(i,j) Y(i,j)];
%                 w(count)=Wx(i)*Wy(j);
%             end
%         end  
%     else
% %         THIS CODE LETS US SELECT THE DESIRED NUMBER OF INTEGRATION POINTS
%         N=nDeg;                %2*nDeg;%fix(sqrt(nDeg));
% %         [X,Y,Wx,Wy]=triquad(N,[-1 -1; 1 -1; -1 1]);
%         [X,Y,Wx,Wy]=triquad(N,[coordTri(1,:); coordTri(2,:); coordTri(3,:) ]);
%         z=zeros(N*N,2);
%         w=zeros(1,N*N);
%         count=0;
%         for i=1:N
%             for j=1:N
%                 count=count+1;
%                 z(count,:)=[X(i,j) Y(i,j)];
%                 w(count)=Wx(i)*Wy(j);
%             end
%         end  
%     end
% 
%     nIP = length(w); %number of integration points
%     nOfGauss = nIP; % for the cubature
% %     Vandermonde matrix
%     V = Vandermonde_LP(nDeg,coord);
%     [L,U,P] = lu(V');
% 
%     shapeFunctions = zeros(nOfNodes,nOfGauss,3);
%     x = z; % (r,s) coordinates
%     [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg);
%     shapeFunctions(:,:,1) = U\(L\(P*p));
%     shapeFunctions(:,:,2) = U\(L\(P*p_xi));
%     shapeFunctions(:,:,3) = U\(L\(P*p_eta));
%     gaussWeights = w';
%     gaussPoints = x;

%% TETRAHEDRA
% N=nDeg^2;       
% vert = [ 1 -1 -1
%         -1 -1 -1
%         -1  1 -1
%         -1 -1  1 ];
% 
% [X,Y,Z,W]=tetraquad(N,vert);
% x =  [ X Y Z ];
% 
% V = Vandermonde_LP(nDeg,coord);
% [L,U,P] = lu(V');
% 
% nOfNodes = (nDeg+1)*(nDeg+2)*(nDeg+3)/6;
% nIP = length(W); %number of integration points
% nOfGauss = nIP; 
% shapeFunctions = zeros(nOfNodes,nOfGauss,4);
% 
% [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_rst(x,nDeg);% [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_xietazeta(x,nDeg);
% shapeFunctions(:,:,1) = U\(L\(P*p));
% shapeFunctions(:,:,2) = U\(L\(P*p_xi));
% shapeFunctions(:,:,3) = U\(L\(P*p_eta));
% shapeFunctions(:,:,4) = U\(L\(P*p_zeta));
% gaussWeights = W;
% gaussPoints = x;
% 
% figure(100)
% plot3(x(:,1),x(:,2),x(:,3),'*r')
% pause()

%% Funcio antiga 3D
% nOfNodes = (nDeg+1)*(nDeg+2)*(nDeg+3)/6;%number of nodes/polynomials
% if nOfNodes~=size(coord,1)
%     error('Error computeShapeFunctionsReferenceElement3D')
% end
% 
% nOfGaussPoints = nDeg^2;
% [z,w] = gaussLegendre(nOfGaussPoints,-1,1);
% nIP = length(w); %number of integration points in each direction
% nOfGauss = length(z);%nIP^3;
% 
% %Vandermonde matrix
% V = Vandermonde_LP(nDeg,coord);
% [L,U,P] = lu(V');
% 
% shapeFunctions = zeros(nOfNodes,nOfGauss,4);
% gaussWeights = zeros(nOfGauss,1);
% gaussPoints = zeros(nOfGauss,3);
% 
% iGauss = 1;
% % Integration over [-1,1]^3
% % tic
% % figure(101)
% for i = 1:nIP
%     for j = 1:nIP
%         for k = 1:nIP
%             x = [z(i),z(j),z(k)]; % (r,s,t) coordinates
% % plot3(z(i),z(j),z(k),'*r')
% % hold on
%             [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_rst(x,nDeg);
% %             [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_xietazeta(x,nDeg);
%             N = U\(L\(P*[p,p_xi,p_eta,p_zeta]));
%             shapeFunctions(:,iGauss,1) = N(:,1)';
%             shapeFunctions(:,iGauss,2) = N(:,2)';
%             shapeFunctions(:,iGauss,3) = N(:,3)';
%             shapeFunctions(:,iGauss,4) = N(:,4)';
%             gaussWeights(iGauss) =(w(i)*w(j)*w(k))*((1-x(2))/2)*((1-x(3))/2)^2;
%             % only for PFEM
%             r = x(1); s = x(2); t = x(3);
%             eta = (1/2)*(s-s*t-1-t);
%             xi = -(1/2)*(r+1)*(eta+t)-1;
%             gaussPoints(iGauss,:) = [xi, eta, t];
%             iGauss = iGauss + 1;
% % hold on
% % plot3(xi,eta,t,'ob')
% % hold on
%         end
%     end
% end

% % pause()
% % toc
% % tic
% % x= zeros(nIP^3,3);
% % for i = 1:nIP
% %     for j = 1:nIP
% %         for k = 1:nIP
% %             x(iGauss,:) = [z(i),z(j),z(k)]; % (r,s,t) coordinates
% %             gaussWeights(iGauss) =(w(i)*w(j)*w(k))*((1-x(2))/2)*((1-x(3))/2)^2;
% %              iGauss = iGauss + 1;
% %         end
% %     end
% % end
% % [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_rst(x,nDeg);
% % % [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_xietazeta(x,nDeg);
% % size( U\(L\(P*p)))
% % size(shapeFunctions(:,:,1))
% % shapeFunctions(:,:,1) = U\(L\(P*p));
% % shapeFunctions(:,:,2) = U\(L\(P*p_xi));
% % shapeFunctions(:,:,3) = U\(L\(P*p_eta));
% % shapeFunctions(:,:,4) = U\(L\(P*p_zeta));
% % % % only for PFEM
% % % r = x(1); s = x(2); t = x(3);
% % % eta = (1/2)*(s-s*t-1-t);
% % % xi = -(1/2)*(r+1)*(eta+t)-1;
% % % gaussPoints(iGauss,:) = [xi, eta, t];
% % % iGauss = iGauss + 1;
% % gaussPoints = '';
% % toc



%% AIXO ERA ABANS LA FUNCIO 2D

% function [shapeFunctions,gaussWeights,gaussPoints]=...
%     computeShapeFunctionsReferenceElement2D(nDeg,coord,coordTri)
% 
% %number of nodes/polynomials
% nOfNodes = (nDeg+1)*(nDeg+2)/2;
% if nOfNodes~=size(coord,1)
%     error('Error computeShapeFunctionsReferenceElement2D')
% end
% % if nDeg <12
% % %     NEW CALL FOR CUBATURES
% %     switch nDeg
% %         case 1
% %             OrderCubature = 5; % 7 pdG
% %         case 2
% %             OrderCubature = 10;
% %         case 3
% %             OrderCubature = 10;
% %         case 4
% %             OrderCubature = 15;
% %         case 5
% %             OrderCubature = 15;
% %         case 6
% %             OrderCubature = 15;
% %         case 7
% %             OrderCubature = 15;
% %         case {8,9,10,11}
% %             OrderCubature = 25;
% %     end
% 
%     if( length(coordTri)<=1 )
%         % THIS CODE DOES NOT LET US CHOOSING A LOT OF INTEGRATION POINTS
% %         [z,w] = GaussLegendreCubature2D(OrderCubature);
% %         w = 2*w; z = 2*z -1; %mapping onto the normal reference triangle
%         N=nDeg^2;                    %4*nDeg*nDeg;%nDeg; %--> sembla bastar amb 2*nDeg
%         [X,Y,Wx,Wy]=triquad(N,[-1 -1; 1 -1; -1 1]);
%         z=zeros(N*N,2);
%         w=zeros(1,N*N);
%         count=0;
%         for i=1:N
%             for j=1:N
%                 count=count+1;
%                 z(count,:)=[X(i,j) Y(i,j)];
%                 w(count)=Wx(i)*Wy(j);
%             end
%         end  
%     else
%         % THIS CODE LETS US SELECT THE DESIRED NUMBER OF INTEGRATION POINTS
%         N=nDeg;                %2*nDeg;%fix(sqrt(nDeg));
% %         [X,Y,Wx,Wy]=triquad(N,[-1 -1; 1 -1; -1 1]);
%         [X,Y,Wx,Wy]=triquad(N,[coordTri(1,:); coordTri(2,:); coordTri(3,:) ]);
%         z=zeros(N*N,2);
%         w=zeros(1,N*N);
%         count=0;
%         for i=1:N
%             for j=1:N
%                 count=count+1;
%                 z(count,:)=[X(i,j) Y(i,j)];
%                 w(count)=Wx(i)*Wy(j);
%             end
%         end  
%     end
% 
%     nIP = length(w); %number of integration points
%     nOfGauss = nIP; % for the cubature
%     %Vandermonde matrix
%     V = Vandermonde_LP(nDeg,coord);
%     [L,U,P] = lu(V');
% 
%     shapeFunctions = zeros(nOfNodes,nOfGauss,3);
%     x = z; % (r,s) coordinates
%     [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg);
%     shapeFunctions(:,:,1) = U\(L\(P*p));
%     shapeFunctions(:,:,2) = U\(L\(P*p_xi));
%     shapeFunctions(:,:,3) = U\(L\(P*p_eta));
%     gaussWeights = w';
%     gaussPoints = x;
% 
% %     figure(10)
% %     plot(x(:,1),x(:,2),'*')
% %     pause()
% 
% % else %--> AIXO NO VAL PER TRIANGLES!!!! CONSIDERES QUADRATS!!!!!!
% %     nOfGaussPoints =( nDeg + 2 )*10;
% %     %OLD CALL
% %     [z,w] = gaussLegendre(nOfGaussPoints,-1,1);
% %     nIP = length(w); %number of integration points in each direction
% %     nOfGauss = ((nIP+1)*nIP)/2;% nIP^2;
% % 
% %     %Vandermonde matrix
% %     V = Vandermonde_LP(nDeg,coord);
% %     [L,U,P] = lu(V');
% % 
% %     shapeFunctions = zeros(nOfNodes,nOfGauss,3);
% %     gaussWeights = zeros(nOfGauss,1);
% % %     gaussPoints = zeros(nOfGauss,2);
% % 
% %     iGauss = 1;
% %     %Integration over [-1,1]^2 --> quadrilateral.. we want tri
% %     x=zeros( nOfGauss ,2);
% %     for i = 1:nIP
% % %         for j = nIP:-1:i % 1:nIP
% % %             x(iGauss,:) = [z(i),z(j)]; % (r,s) coordinates
% % %             gaussWeights(iGauss) =(w(i)*w(j))*(1-x(iGauss,2))/2;
% % %             iGauss = iGauss + 1;
% % 
% %         for j = 1:nIP
% %             if( (z(i)+z(j)) <= 0 )
% %                 x(iGauss,:) = [z(i),z(j)]; % (r,s) coordinates
% %                 gaussWeights(iGauss) =(w(i)*w(j))*(1-x(iGauss,2))/2;
% %                 iGauss = iGauss + 1;
% %             end
% %             
% %         end
% %     end
% %     
% % %     figure(10)
% % %     plot(x(:,1),x(:,2),'*')
% % %     pause()
% %     
% %     [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg);
% %     shapeFunctions(:,:,1) = U\(L\(P*p));
% %     shapeFunctions(:,:,2) = U\(L\(P*p_xi));
% %     shapeFunctions(:,:,3) = U\(L\(P*p_eta));
% %     % only for PFEM
% %     r = x(:,1); s = x(:,2);
% %     xi = (1+r).*(1-s)/2-1;
% %     gaussPoints = [xi, s];
% %     
% % 
% % % % Initial code:    
% % %     for i = 1:nIP
% % %         for j = 1:nIP
% % %             x = [z(i),z(j)]; % (r,s) coordinates
% % %             [p,p_xi,p_eta] = orthopoly2D_deriv_rst(x,nDeg);
% % %             N = U\(L\(P*[p,p_xi,p_eta]));
% % %             shapeFunctions(:,iGauss,1) = N(:,1)';
% % %             shapeFunctions(:,iGauss,2) = N(:,2)';
% % %             shapeFunctions(:,iGauss,3) = N(:,3)';
% % %             gaussWeights(iGauss) =(w(i)*w(j))*(1-x(2))/2;
% % %             % only for PFEM
% % %             r = x(1); s = x(2);
% % %             xi = (1+r)*(1-s)/2-1;
% % %             gaussPoints(iGauss,:) = [xi, s];
% % %             iGauss = iGauss + 1;
% % %         end
% % %     end
% % end



%