function []=function_createBlocks(H_build,W_build,factor_gap_from_H)

%W_gap,H_build,W_build

if(nargin==0)
    H_build = 30;
    W_build = [H_build H_build];
    factor_gap_from_H = [0 1 2 3];
end

factor_gap_from_H = [0 factor_gap_from_H];

numBlocks = length(factor_gap_from_H)

%fileName = '3BlocksSynthetic';
fileName = ['blocks_N' int2str(numBlocks)];



X = zeros(2,0);
T = zeros(2,0);


Xcorner = [0.0 ;0.0];

for ifactor = 1:numBlocks
    
    factorH = factor_gap_from_H(ifactor);
    
    if(ifactor>1)
        fileName = [fileName '_' int2str(factorH)]
    end
    
    W_gap = H_build*factorH;
    
    x0 = Xcorner + [W_gap;0];
    
    [Xblock,Tblock] = createBlock(x0,W_build);
    
    T = [T (Tblock+size(X,2))];
    X = [X Xblock];
    
    Xcorner = Xblock(:,2);
    
end

ycm = sum(X(2,:))/size(X,2);
X(2,:) = X(2,:)-ycm;

figure(1)
plot(X(1,:),X(2,:),'*')
axis equal

exportBlocks(X,T,fileName)

exportTranslationPoint(fileName)

exportField(H_build,[fileName '_lidar'])

exportField(0.0,[fileName '_dem'])


end

function [Xblock,Tblock]=createBlock(x0,W)

    if(length(W)==1)
        Xblock = [x0,x0+[W;0],x0+[W;W],x0+[0;W]];
    else
        Xblock = [x0 , x0+[W(1);0] , x0+[W(1);W(2)] , x0+[0;W(2)]];
    end
    Tblock = [ 1 2 3 4
               2 3 4 1];
    
end

function []=exportBlocks(X,T,fileName)
    
    numNodes = size(X,2);
    numElems = size(T,2);

    stream = fopen(['./input/' fileName '.msh'],'w');
    % write header
    fprintf(stream,['MESH dimension 3 ElemType Linear Nnode 2\n']);
    % write node coordinates
    fprintf(stream,['Coordinates\n']);
    fprintf(stream,'%d %f %f 0\n', [1:numNodes; X]);
    fprintf(stream,['End Coordinates\n\n']);
    % write elements
    fprintf(stream,['Elements\n']);
    fprintf(stream,'%d %d %d\n', [ 1:numElems; T]);
    fprintf(stream,['End Elements\n']);
    %
    fclose(stream);

end

function []=exportTranslationPoint(fileName)
    
    stream = fopen(['./input/' fileName '_translationPoint.txt'],'w');
    fprintf(stream,['0.0 0.0\n']);
    fclose(stream);

end

function []=exportField(H,fileName)
    
    hx = 1e10;
    hy = hx;
    x0 = [-hx  -hy];
    nx = 2;
    ny = 2;
    
    field.structured = 1;
    field.z = H*ones(nx,ny);
    field.hx = hx;
    field.hy = hy;
    field.x0 = x0;
    field.nx = nx;
    field.ny = ny;

    save(['./input_fields/' fileName '.mat'],'field');
    

end

% MESH dimension 3 ElemType Linear Nnode 2
% Coordinates
%     1     -900.911185     -274.876506               0
%     2      -900.82219     -254.868701               0
%     3     -962.218576      -335.12693               0
%     4     -982.353381     -334.780938               0
%     5     -961.312615     -193.221315               0
%     6     -981.506419     -193.291319               0
%     7     -1042.75381     -273.127552               0
%     8     -1042.67181     -253.642743               0
% End Coordinates
% 
% Elements
% 1 3 4
% 2 2 1
% 3 1 3
% 4 4 7
% 5 7 8
% 6 8 6
% 7 6 5
% 8 5 2
% End Elements
