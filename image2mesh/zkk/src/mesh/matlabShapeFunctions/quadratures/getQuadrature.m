function [x,w] = getQuadrature(element,N)

    switch element.type
        
        case 'tri'
            [x,w] = getQuadrature_tri(element,N);
            
        case 'quad'
            [x,w] = getQuadrature_quad(element,N);

        case 'tet'
            [x,w] = getQuadrature_tet(element,N);
            
        case 'hex'
            [x,w] = getQuadrature_hex(element,N);

        otherwise 
            error('this function is not yet implemented for this element')
            
    end

end

function [x,w] = getQuadrature_hex(element,N)

    [x,w] = lglnodesHex(N);
    
end

function [x,w] = getQuadrature_quad(element,N)

        [x,w] = lglnodesQuad(N);

        % to send to another quad I must use what I will program now.
%         quad1 = [ [-1;-1] [1;-1] [1;1] [-1;1] ];
%         quad2 = coordElem(permutes(iPerm,:),:)';
%         
%         [x] = sendQuad1PointsToQuad2(x',quad1,quad2);
%         x = x';
        

end

function [x,w] = getQuadrature_tri(element,N)
   
vert = element.vert;

if(isfield(element,'quadrature'))
    switch element.quadrature

        case 'gaussLobatto'
%             fprintf('Gauss-lobatto quadrature\n')
             [x,w]=lglnodes_tri(vert',N);
%              w = ones(length(w),1);

        case  'gauss'
%             fprintf('Gauss quadrature\n')
%             disp('integracio triquad')
            %permutes = [1 2 3]; % just one quadrature
            % permutes = [ 1 2 3 ; 3 1 2 ; 1 3 2 ]; % symmetrized on vertices, not on edges    
            permutes = perms([1 2 3 ]); % fully symmetrized     
            numPerms = size(permutes,1);

            nOfGauss = N^2; 
            numOfGaussSym   = nOfGauss*numPerms;

             x = zeros(numOfGaussSym,2);
             w = zeros(numOfGaussSym,1);
            
            count=0;
            for iPerm = 1:numPerms
                [X,Y,Wx,Wy]=triquad(N,vert(permutes(iPerm,:),:));
                %[X,Y,Wx,Wy]=triquad_lobatto(N,coordElem(permutes(iPerm,:),:));%no correcte per lobatto
                for i=1:N
                    for j=1:N
                        count=count+1;
                        x(count,:)=[X(i,j) Y(i,j)];
                        w(count)=Wx(i)*Wy(j);
                    end
                end  
                %[x,w]=triquadReduced_lobatto(N,coordElem(permutes(iPerm,:),:));
            end
        case  'gaussStandard'
%             fprintf('Gauss quadrature\n')
%             disp('integracio triquad')
            permutes = [1 2 3]; % just one quadrature
            % permutes = [ 1 2 3 ; 3 1 2 ; 1 3 2 ]; % symmetrized on vertices, not on edges    
            %permutes = perms([1 2 3 ]); % fully symmetrized     
            numPerms = size(permutes,1);

            nOfGauss = N^2; 
            numOfGaussSym   = nOfGauss*numPerms;

             x = zeros(numOfGaussSym,2);
             w = zeros(numOfGaussSym,1);
            
            count=0;
            for iPerm = 1:numPerms
                [X,Y,Wx,Wy]=triquad(N,vert(permutes(iPerm,:),:));
                %[X,Y,Wx,Wy]=triquad_lobatto(N,coordElem(permutes(iPerm,:),:));%no correcte per lobatto
                for i=1:N
                    for j=1:N
                        count=count+1;
                        x(count,:)=[X(i,j) Y(i,j)];
                        w(count)=Wx(i)*Wy(j);
                    end
                end  
                %[x,w]=triquadReduced_lobatto(N,coordElem(permutes(iPerm,:),:));
            end
    end
end
end

function [x,w] = getQuadrature_tet(element,N)

    vert = element.vert;

    if(isfield(element,'quadrature'))
        switch element.quadrature                
            
            case 'gaussLobatto'
%                  fprintf('Gauss-lobatto quadrature\n')
                 [x,w]=lglnodes_tet(vert',N);
%                  w = ones(length(w),1);
%                   w = 1./w;
            
            case 'gaussSym'
%                 fprintf('Gauss sym using subhex')
                [x,w]=gaussSym_tet(vert',N);

            case  'gauss'
%                 fprintf('Gauss integration points \n')
                
                permutes = [1 2 3 4]; % just one quadrature
                % permutes = [ 1 2 3 4 ; 3 1 2 4 ; 1 3 2 4 ; 1 4 2 3 ]; % symmetrized on vertices, not on edges    
                % permutes = perms([1 2 3 4]); % fully symmetrized     
                numPerms = size(permutes,1);
                
                nOfGauss = N^3; 
                numOfGaussSym   = nOfGauss*numPerms;
                
                 x = zeros(numOfGaussSym,3);
                 w = zeros(numOfGaussSym,1);
                
                for iPerm = 1:numPerms
                    iGauss_ini = 1 + nOfGauss*(iPerm-1);
                    iGauss_end = iGauss_ini + nOfGauss - 1;
                    
                    [X,Y,Z,W]=tetraquad(N,vert(permutes(iPerm,:),:));
                    x(iGauss_ini:iGauss_end,:) = [ X Y Z ];
                    w(iGauss_ini:iGauss_end) = W/numPerms;
                end
                    
            case 'nodes'
%                 fprintf('Node as integration points \n')
                x = giveReferencePoints_tet(element);
                w = 1/size(x,1)*ones(size(x,1),1);
                
            case 'nodesHigherOrder'
%                 fprintf('Node as integration points \n')
                dim=3;type='tet'; quadrature ='';
                distribution = 'fekete'; % 'fekete', 'equispaced', 'lineal'
                higherOrder = element.order*2;
                elementHigherOrder=defineElement(dim,dim,type,higherOrder,distribution,quadrature);
                x = giveReferencePoints_tet(elementHigherOrder);
                w = 1/size(x,1)*ones(size(x,1),1);    

            otherwise
                error('quadrature not implemented')
                
        end
    else
        error('specify quadrature for tets');
    end


end