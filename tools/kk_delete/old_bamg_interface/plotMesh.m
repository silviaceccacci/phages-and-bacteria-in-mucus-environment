function plotMesh(X,T,colouring,Thigh,feketeDistribution)

    if(nargin<4 || Thigh(1,1)==0)

        if(length(colouring)==1 && colouring(1)==0)
            if(size(X,2)<3)
%                 trisurf(T, X(:,1), X(:,2), zeros(size(X(:,1))) )
                  triplot(T, X(:,1), X(:,2))
            else
                trisurf(T, X(:,1), X(:,2), X(:,3))
            end
            
        else
            
            if(size(X,2)<3)
                trisurf(T, X(:,1), X(:,2), zeros(size(X(:,1))) , colouring)
            else
                trisurf(T, X(:,1), X(:,2), X(:,3), colouring)
            end            
            colorbar
            caxis([0 1])
        end

        
    else
        
        if(length(colouring)==1 && colouring(1)==0)

            if(size(X,2)<3)
                trisurf(T, X(:,1), X(:,2), zeros(size(X(:,1))) ,'EdgeColor','none')
            else
                trisurf(T, X(:,1), X(:,2), X(:,3),'EdgeColor','none')
            end
            
        else

            
            if(size(X,2)<3)
                trisurf(T, X(:,1), X(:,2), zeros(size(X(:,1))) , colouring,'EdgeColor','none')
            else
                trisurf(T, X(:,1), X(:,2), X(:,3), colouring,'EdgeColor','none')
            end
            colorbar
            caxis([0 1])

        end

        plotMesh_HO(X,Thigh,feketeDistribution);  
        
    end

%     list = 1:size(X,1);
%     fontSize = 15;
%     for inode = list
%        text(X(inode,1),X(inode,2),int2str(inode),'FontSize',fontSize)
%     end
    
    
    set(gca,'Fontsize',16); 
    colormap(jet)
    lighting flat

    if(size(X,2) < 3)
        view(2)
    else
        numZ0=length(find(X(:,3)==0));
        if(numZ0==size(X,1))
            view(2)
        end  
    end

    axis equal
    
end



% if(dimension==3)
%     if(length(colouring)==1 && colouring(1)==0)
%         trisurf(T, X(:,1), X(:,2), X(:,3))
%     else
%         trisurf(T, X(:,1), X(:,2), X(:,3), colouring)
%         colorbar
%         caxis([0 1])
%     end
% 
% elseif(dimension==2)
%     if(length(colouring)==1 && colouring(1)==0)
%         triplot( T, X(:,1), X(:,2) )
%     else
% %         triplot( T, X(:,1), X(:,2) )%, colouring )
%         trisurf(T, X(:,1), X(:,2), X(:,3), colouring)
%         colorbar
%         caxis([0 1])
%         axis equal;
%     end   
% else
%     disp('You have to fix a plot dimension: 2D or 3D plot')
% end
% 
% 
% set(gca,'Fontsize',16); 
% colormap(jet)
% lighting flat
