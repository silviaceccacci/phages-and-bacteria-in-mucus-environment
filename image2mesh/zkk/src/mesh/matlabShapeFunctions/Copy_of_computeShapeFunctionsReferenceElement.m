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
% nOfGaussPoints: n� of gauss points of the 1D quadrature
% elementType (optional): 0 for quadrilateral, 1 for triangle. If it isn't
%                         given only triangle or 1D elements are considered
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

[z,w] = gaussLegendre(nOfGaussPoints,-1,1);
nOfGauss = length(w);

%Vandermonde matrix
V = Vandermonde_LP(nDeg,coord);
[L,U,P] = lu(V');

shapeFunctions = zeros(nOfNodes,nOfGauss,2);
gaussWeights = w';
% gaussPoints = zeros(nOfGauss,1);

%Integration over [-1,1]
x = z' ;
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
nOfNodes = giveNumNodesElementFromOrder(nDeg,element.type);
if nOfNodes~=size(coord,1)
    error('Error computeShapeFunctionsReferenceElement2D')
end
switch element.type
    case 'tri'
        if( length(coordElem)<=1 || isempty( find( coordElem~= [-1 -1; 1 -1; -1 1])  ) )
            N = nDeg+1;%     N = 2*nDeg +1;
            coordElem = [-1 -1; 1 -1; -1 1];
        else
            N = ceil(sqrt(nDeg));%     N= 2*ceil(sqrt(nDeg))+1;
        end
% % %-----------    
        permutes = perms([1 2 3]);
%         permutes = [1 2 3 ; 2 3 1 ; 3 1 2 ]; % fully symmetrized (for tris)
%         permutes = [1 2 3]; %disp('not symmetric integration')
% % %-----------
%         permutes = 0;
%         switch nDeg
%             case 1
%                 OrderCubature = 5; % 7 pdG
%             case {2,3}
%                 OrderCubature = 10;
%             case {4,5,6,7}
%                 OrderCubature = 15;
%             case {8,9,10,11}
%                 OrderCubature = 25;
%         end
%         [z,w] = GaussLegendreCubature2D(OrderCubature);
%         w = 2*w; x = 2*z -1; %mapping onto the normal reference triangle
%         N =0;
%         disp('integracio quadratura paper')
% % %-----------
    case 'quad'
        if( length(coordElem)<=1 || isempty( find( coordElem~= [-1 -1; 1 -1; 1 1 ;-1 1])  ) )
            N = nDeg+2;
            coordElem = [-1 -1; 1 -1; 1 1 ;-1 1];
        else
            N = ceil(sqrt(nDeg))+2;
        end
        permutes = [1 2 3 4]; %symmetric (since it is default symmetric)
end
numPerms = size(permutes,1);

V = Vandermonde_LP(nDeg,coord,element);
[L,U,P] = lu(V');

if(N>0)
    nOfGauss = N^2;
else
    nOfGauss = length(w);
end
% [x,w] = triquadReduced_lobatto(N);
% nOfGauss = length(w);

numOfGaussSym   = nOfGauss*numPerms;
shapeFunctions  = zeros(nOfNodes,numOfGaussSym,3);
gaussWeights    = zeros(numOfGaussSym,1);
gaussPoints     = zeros(numOfGaussSym,2);

for iPerm = 1:numPerms
    iGauss_ini = 1 + nOfGauss*(iPerm-1);
    iGauss_end = iGauss_ini + nOfGauss - 1;
    
    if(strcmp(element.type,'tri'))
        if(length(permutes)>1)
            disp('integracio triquad')
            [X,Y,Wx,Wy]=triquad(N,coordElem(permutes(iPerm,:),:));
% %     %         [X,Y,Wx,Wy]=triquad_lobatto(N,coordElem(permutes(iPerm,:),:)); %    no correcte per lobatto
            x=zeros(nOfGauss,2);
            w=zeros(1,nOfGauss);
            count=0;
            for i=1:N
                for j=1:N
                    count=count+1;
                    x(count,:)=[X(i,j) Y(i,j)];
                    w(count)=Wx(i)*Wy(j);
                end
            end  
%             disp('integarcio triquad lobatto')
%             [x,w] = triquadReduced_lobatto(N,coordElem(permutes(iPerm,:),:));
        end
        
        [p,p_xi,p_eta] = orthopoly2D_deriv_xieta(x,nDeg,element); % AQUI ES XIETA I A 3D RST!?
        
    elseif(strcmp(element.type,'quad'))
        
        [x,w] = lglnodesQuad(N);

        quad1 = [ [-1;-1] [1;-1] [1;1] [-1;1] ];
        quad2 = coordElem(permutes(iPerm,:),:)';
        
        [x] = sendQuad1PointsToQuad2(x',quad1,quad2);
        x = x';
        
     % % Els punts de Gauss i els de Fekete coincideixen!   
%     [coordGauss,kk] = lglnodesQuad(element.order+1);    
%     [ coordFekete ] = generateFeketeElementFromVertices_quad(element);
%     figure(9)
%     clf
%     plot(coordGauss(:,1),coordGauss(:,2),'b*','markersize',15)
%     hold on
%     plot(coordFekete(:,1),coordFekete(:,2),'r*')
    
    % % Pintem els punts de gauss de la regio dintegracio
%     figure(10)
%     clf
%     plot(quad1(1,:),quad1(2,:),'o','markersize',10)
%     hold on
%     plot(quad2(1,:),quad2(2,:),'ro','markersize',5)
%     hold on
%     plot(x(:,1),x(:,2),'b*','markersize',5)
%     pause()

        
        [p,p_xi,p_eta] = orthopoly2D_deriv_quads(x,nDeg,element);
        
    else
        
        error('Error computeShapeFunctionsReferenceElement2D in element type')
        
    end     
    shapeFunctions(:,iGauss_ini:iGauss_end,1) = U\(L\(P*p));
    shapeFunctions(:,iGauss_ini:iGauss_end,2) = U\(L\(P*p_xi));
    shapeFunctions(:,iGauss_ini:iGauss_end,3) = U\(L\(P*p_eta));
    gaussWeights(iGauss_ini:iGauss_end) = w/numPerms;
    gaussPoints(iGauss_ini:iGauss_end,:) = x;
end
% plot(gaussPoints(:,1),gaussPoints(:,2),'*')
% pause()
% shapeFunctions
% pause()

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement3D(nDeg,coord,options,element)
%% SYMMETRIZED QUADRATURE 
if(isfield(options,'vert')==false || max( isempty(options.vert) ) ||  (options.vert(1,1) == 0) || ...
        isempty( find( options.vert~= [ 1 -1 -1;-1 -1 -1;-1  1 -1;-1 -1  1 ])  ) )
    vert = [ 1 -1 -1
            -1 -1 -1
            -1  1 -1
            -1 -1  1 ];
%     N = 2*nDeg+1; %     N = nDeg;
    N = nDeg+1;%nDeg+1;
else
    vert = options.vert;
    N = 2*ceil(sqrt(nDeg)) +1;
end

permutes = [1 2 3 4]; % just one quadrature
% permutes = [ 1 2 3 4 ; 3 1 2 4 ; 1 3 2 4 ; 1 4 2 3 ]; % symmetrized on vertices, not on edges    
% permutes = perms([1 2 3 4]); % fully symmetrized     

numPerms = size(permutes,1);

nOfNodes = (nDeg+1)*(nDeg+2)*(nDeg+3)/6;
nOfGauss = N^3; 
% nOfGauss = nOfNodes;
numOfGaussSym   = nOfGauss*numPerms;
shapeFunctions  = zeros(nOfNodes,numOfGaussSym,4);
gaussWeights    = zeros(numOfGaussSym,1);
gaussPoints     = zeros(numOfGaussSym,3);

V = Vandermonde_LP(nDeg,coord,element);
[L,U,P] = lu(V');

for iPerm = 1:numPerms
    iGauss_ini = 1 + nOfGauss*(iPerm-1);
    iGauss_end = iGauss_ini + nOfGauss - 1;

    if(isfield(element,'quadrature'))
        switch element.quadrature
            case  'gauss'
                fprintf('\n gauss integration points \n')
                [X,Y,Z,W]=tetraquad(N,vert(permutes(iPerm,:),:));
                x =  [ X Y Z ];
            case 'nodes'
                fprintf('\n node as integration points \n')
                x = giveReferencePoints_tet(element);
                W = 1/size(x,1)*ones(size(x,1),1);
                nOfGauss = element.numNod;
                iGauss_end = iGauss_ini + nOfGauss - 1;
            case 'nodesHigherOrder'
                fprintf('\n node as integration points \n')
                dim=3;type='tet'; quadrature ='';
                distribution = 'fekete'; % 'fekete', 'equispaced', 'lineal'
                higherOrder = element.order*2;
                elementHigherOrder=defineElement(dim,dim,type,higherOrder,distribution,quadrature);
                nOfGauss = elementHigherOrder.numNod;
                iGauss_end = iGauss_ini + nOfGauss - 1;
                x = giveReferencePoints_tet(elementHigherOrder);
                W = 1/size(x,1)*ones(size(x,1),1);    
            case 'gaussLobatto'
                 [x,W]=lglnodes_tet(vert',N);
            otherwise
                error('quadrature not implemented')
        end
    else
        error('specify quadrature for tets');
    end

    [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_rst(x,nDeg);% [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_xietazeta(x,nDeg);

    shapeFunctions(:,iGauss_ini:iGauss_end,1) = U\(L\(P*p));
    shapeFunctions(:,iGauss_ini:iGauss_end,2) = U\(L\(P*p_xi));
    shapeFunctions(:,iGauss_ini:iGauss_end,3) = U\(L\(P*p_eta));
    shapeFunctions(:,iGauss_ini:iGauss_end,4) = U\(L\(P*p_zeta));
    gaussWeights(iGauss_ini:iGauss_end) = W/numPerms;
    gaussPoints(iGauss_ini:iGauss_end,:) = x;
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