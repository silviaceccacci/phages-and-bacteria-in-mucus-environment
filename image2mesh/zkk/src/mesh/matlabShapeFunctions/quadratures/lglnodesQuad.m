function [x,w,X,Y,Wx,Wy] = lglnodesQuad(N)
%%%%%%%%%%%%%%%%%%%%%%
% Created by AGP
% This function produces N^2 Gauss-Lobatto 
% nodes on a reference quadrilateral element
% [-1,1] x [-1,1]
%%%%%%%%%%%%%%%%%%%%%%

    [x,w]  = lglnodes(N-1);

    X        = repmat(x,1,length(x));
    Wx     = repmat(w,1,length(x));
    Y        = X';
    Wy     = Wx';

    
    numIntPoints = length(x)^2;
    x  = zeros(numIntPoints,2);
    w = zeros(numIntPoints,1);
    count = 0;
    for i=1:size(X,1)
        for j=1:size(X,2)
            count=count+1;
            x(count,:)=[X(i,j) Y(i,j)];
            w(count)=Wx(i,j)*Wy(i,j);
        end
    end  
    
end