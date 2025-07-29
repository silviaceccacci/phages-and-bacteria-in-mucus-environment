function [xlgl,wlgl]=gaussSym_tet(tetVert,N)

    % define quadrature in an hex
    [xlgl_hex,weights_hex] = gaussHex(N);

    % define shape functions on an hex of order 1 on xlg_hex
    [shapeFunctions] = defineShapeFunctionsHexOrder1(xlgl_hex);

    % define subhexes
    [Xsubhex,Tsubhex] = defineSubHexesFromTet(tetVert);
    numSubHexes = size(Tsubhex,1);

    % map quadrature to the four subhexes
    numIntPointsHex = size(xlgl_hex,1);
    xlgl = zeros(numIntPointsHex*numSubHexes,size(xlgl_hex,2));
    wlgl = zeros(numIntPointsHex*numSubHexes,1);
    iNew = 1;
    for iSubHex = 1:numSubHexes
        iNewFinal = iNew + numIntPointsHex-1;
        
        xISubHex = Xsubhex(:,Tsubhex(iSubHex,:));  
        xlgl(iNew:iNewFinal,:) = (xISubHex*shapeFunctions(:,:,1))';

        J = zeros(size(shapeFunctions,2),size(xISubHex,1),3);
        J(:,:,1) =  ( xISubHex*shapeFunctions(:,:,2) )';
        J(:,:,2) =  ( xISubHex*shapeFunctions(:,:,3) )';
        J(:,:,3) =  ( xISubHex*shapeFunctions(:,:,4) )';
        detJ = determinantNDim(J);
        wlgl(iNew:iNewFinal) = abs(detJ).*weights_hex;
% figure(2)
% clf;
% % plot3(xlgl_hex(:,1),xlgl_hex(:,2),xlgl_hex(:,3),'*')
% hold on
% plot3(Xsubhex(1,:),Xsubhex(2,:),Xsubhex(3,:),'*')
% plot3(xISubHex(1,:),xISubHex(2,:),xISubHex(3,:),'ro')
% axis equal
% sum(weights_hex)
% sum(wlgl(iNew:iNewFinal))
% pause()
% %         xlgl(iNew:iNewFinal,:)
%         figure(24)
%         plot3(xlgl_hex(:,1),xlgl_hex(:,2),xlgl_hex(:,3),'*')
%         figure(25);  clf; axis equal; hold on
%         plot3(Xsubhex(1,:),Xsubhex(2,:),Xsubhex(3,:),'*')
%         plot3(xISubHex(1,:),xISubHex(2,:),xISubHex(3,:),'o')
%         plot3(xlgl(iNew:iNewFinal,1),xlgl(iNew:iNewFinal,2),xlgl(iNew:iNewFinal,3),'ro','MarkerSize',10)
%         pause()

        iNew = iNewFinal + 1;
    end 
    
end

function [shapeFunctions] = defineShapeFunctionsHexOrder1(xlgl_hex)

        dim=3;type='hex'; quadrature =''; distribution = 'linear'; % 'fekete', 'equispaced', 'lineal'
        elementHex=defineElement(dim,dim,type,1,distribution,quadrature);

        coordRefHex = [ 1  -1  -1
                                1   1  -1
                               -1   1   -1
                               -1  -1  -1
                                1  -1   1
                                1   1    1
                                -1  1    1
                                -1  -1   1];
%                 coordRefHex = (coordRefHex +1)/2; % in [0,1]^3

        V = Vandermonde_LP(elementHex.order,coordRefHex,elementHex);
        [L,U,P] = lu(V');

        [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_hex(xlgl_hex,elementHex);
        
        numIntPointsHex = size(xlgl_hex,1);
        shapeFunctions  = zeros(elementHex.numNod,numIntPointsHex,4);
        shapeFunctions(:,:,1) = U\(L\(P*p));
        shapeFunctions(:,:,2) = U\(L\(P*p_xi));
        shapeFunctions(:,:,3) = U\(L\(P*p_eta));
        shapeFunctions(:,:,4) = U\(L\(P*p_zeta));

end


function [Xsubhex,Tsubhex] = defineSubHexesFromTet(tetVert)

        Tsubhex = [ 1  2  3  4  5  6  7  8
                        9 10 3  2 11 12 7 6 
                        13 4 3 10 14 8 7 12
                        15 5 8 14 11 6 7 12];
        numSubNod = 15;
        Xbar = [   1       0       0     0       %1
                      1/2    1/2     0     0       %2
                      1/3    1/3    1/3   0       %3
                       1/2    0    1/2   0       %4
                       1/2    0       0     1/2    %5
                       1/3   1/3     0     1/3    %6
                       1/4   1/4    1/4  1/4     %7
                       1/3    0      1/3   1/3    %8
                       0       1       0      0       %9
                       0       1/2    1/2   0       %10
                       0       1/2    0      1/2    %11
                       0       1/3    1/3   1/3    %12
                       0       0        1      0       %13
                       0       0        1/2   1/2    %14
                       0       0         0      1  ];  %15
        Xsubhex = zeros(3,numSubNod);
        for iSubNod = 1:numSubNod
            Xsubhex(:,iSubNod) = tetVert*Xbar(iSubNod,:)';
        end

%         figure(50)
%         plot3(Xsubhex(1,:),Xsubhex(2,:),Xsubhex(3,:),'*')
%        axis equal
% for i=1:4
%         figure(55)
%         clf;
%         plot3(Xsubhex(1,:),Xsubhex(2,:),Xsubhex(3,:),'*')
%         hold on
%         plot3(Xsubhex(1,Tsubhex(i,:)),Xsubhex(2,Tsubhex(i,:)),Xsubhex(3,Tsubhex(i,:)),'ro')
%         axis equal
%         pause()
% end


end
