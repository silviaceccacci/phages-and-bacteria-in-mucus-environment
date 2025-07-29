function [boolea]=isBoundary(Telem,boundaryNodes)

for i=1:3%length(Telem)
    
    if( isempty(find(Telem(i))==boundaryNodes) == false )
    
        boolea=true;
        
        return
        
    end
    
end

boolea=false;
        