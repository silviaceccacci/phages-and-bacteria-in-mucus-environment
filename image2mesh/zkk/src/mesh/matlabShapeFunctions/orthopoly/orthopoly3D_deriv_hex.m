function [p,p_xi,p_eta,p_zeta] = orthopoly3D_deriv_hex(x,element)

    nDeg = element.order;

    [p_x,dp_x]  = orthopoly1D_deriv(x(:,1), nDeg); 
    [p_y,dp_y]  = orthopoly1D_deriv(x(:,2), nDeg);
    [p_z,dp_z]  = orthopoly1D_deriv(x(:,3), nDeg);

    numNod = element.numNod;
    p = zeros(numNod, size(p_x,2)); 
    p_xi = p; p_eta = p; p_zeta = p;
    count = 0;
    for ix = 1:size(p_x,1) %nDeg+1
        for iy = 1:size(p_y,1) %nDeg+1
            for iz = 1:size(p_z,1) %nDeg+1
                count = count+1;
                p(count,:)        =  p_x(ix,:).*   p_y(iy,:) .*  p_z(iz,:) ;
                p_xi(count,:)   =  dp_x(ix,:).*  p_y(iy,:) .*  p_z(iz,:) ;
                p_eta(count,:) =   p_x(ix,:).* dp_y(iy,:) .*  p_z(iz,:) ;
                p_zeta(count,:) =  p_x(ix,:).*  p_y(iy,:) .* dp_z(iz,:) ;
            end
        end
    end  
    
end

% function [p,p_xi,p_eta] = orthopoly2D_deriv_hex(x,element)
% %% still on development
% %%%%%%%%%%%%%%%%%%%%%%
% % INPUT
% % size(x)  -->gauss x dim
% %%%%%%%%%%%%%%%%%%%%%%
% 
%     nDeg = element.order;
%     n = element.numNod;
% 
%     xFactor = 1e6;
%     xMod = round( x * xFactor ); 
%     coord_x = unique(xMod(:,1)) ./xFactor;
%     coord_y = unique(xMod(:,2)) ./xFactor;
%     coord_z = unique(xMod(:,3)) ./xFactor;
%     
%     [p_x,dp_x]  = orthopoly1D_deriv(coord_x, nDeg); % returns matrix (sqrt(n),size(coord_x,1))
%     [p_y,dp_y]  = orthopoly1D_deriv(coord_y, nDeg); % returns matrix (sqrt(n),size(coord_x,1))
%     [p_z,dp_z]  = orthopoly1D_deriv(coord_z, nDeg); % returns matrix (sqrt(n),size(coord_x,1))
%    
%     p        = zeros(n,size(x,1)); % size(p) --> nodes x gauss
%     p_xi   = zeros(n,size(x,1)); % size(p_xi) --> nodes x gauss
%     p_eta = zeros(n,size(x,1)); % size(p_eta) --> nodes x gauss
%     p_zeta = zeros(n,size(x,1)); % size(p_eta) --> nodes x gauss
%     count = 0;
%     for ix = 1:size(p_x,1) %nDeg+1
%         for iy = 1:size(p_y,1) %nDeg+1
%             for iz = 1:size(p_z,1) %nDeg+1
%                 count = count+1;
%                 p(count,:)        = reshape((   p_x(ix,:)'*  p_y(iy,:) ), 1, size(p,2) );
%                 p_xi(count,:)   = reshape(( dp_x(ix,:)'*  p_y(iy,:) ), 1, size(p,2) );
%                 p_eta(count,:) = reshape((   p_x(ix,:)'*dp_y(iy,:) ), 1, size(p,2) );
%                 p_zeta(count,:) = reshape((   p_x(ix,:)'*dp_y(iy,:) ), 1, size(p,2) );
%             end
%         end
%     end
%     
%     %% Reorder node location     
%     coordRef = giveReferencePoints_quad(element);
%     [coordRefSorted , index0] = sortrows([coordRef(:,1)  coordRef(:,2)]);
%     index = zeros(size(index0));
%     for iaux = 1:length(index0)
%         index(iaux) = find(iaux == index0);
%     end
%     p        = p(index,:);
%     p_xi   = p_xi(index,:);
%     p_eta = p_eta(index,:);
%     
%     %% Reorder gauss points locations to be in quadrilateral order
%     % so that it can be added with the corresponding gauss weight
%     % no se si aixo esta be aixi o si iy i ix no han danar girats
%     x2 = zeros(size(x));
%     count = 0;
%     for iy = 1:size(coord_y,1) % sqrt(gauss)
%         for ix = 1:size(coord_x,1) % sqrt(gauss)
%             count = count+1;
%             x2(count,:) = [coord_x(ix) coord_y(iy)];
%         end
%     end
%     
%     listOrder = zeros(size(x,1),1);
%     for aux = 1:size(x,1)
%          ind1 = find(abs(x2(:,1)-x(aux,1))<1/xFactor);
%          listOrder(aux) = ind1( find(abs(x2(ind1,2)-x(aux,2))<1/xFactor) );
%     end
%     p        = p(:,listOrder);
%     p_xi   = p_xi(:,listOrder);
%     p_eta = p_eta(:,listOrder);
%     
% end