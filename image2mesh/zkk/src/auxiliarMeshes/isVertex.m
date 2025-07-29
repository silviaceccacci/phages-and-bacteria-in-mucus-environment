function [ isVertex_var ] = isVertex( localNum , Telem)

%% JUST FOR HIGH ORDER TRIANGLES

if(nargin == 2) % if the input is the global numbering of the node
    localNum = find(Telem == localNum);
end



isVertex_var = true;

if( localNum > 3 )
    
    isVertex_var = false;
    
end