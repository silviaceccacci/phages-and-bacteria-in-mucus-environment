function [x,w] = lglnodesHex(N)
%%%%%%%%%%%%%%%%%%%%%%
% Created by AGP
% This function produces N^2 Gauss-Lobatto 
% nodes on a reference hexahedral element
% [    1 -1 -1
%     -1 -1 -1
%     -1  1 -1
%     -1 -1  1 ];
%%%%%%%%%%%%%%%%%%%%%%

    [x1,w1]  = lglnodes(N-1);
    numIntPoints1 = length(x1);
    
    numIntPoints = numIntPoints1^3;
    x  = zeros(numIntPoints,3);
    w = zeros(numIntPoints,1);
    count = 0;
    for i=1:size(x1,1)
        for j=1:size(x1,1)
            for k=1:size(x1,1)
                count=count+1;
                x(count,:)=[x1(i) x1(j) x1(k)];
                w(count)=w1(i)*w1(j)*w1(k);
            end
%             countOld = count;
%             count=count+numIntPoints1;
%             x(countOld:count,:)=[x1(i)*ones(numIntPoints1,1) x1(j)*ones(numIntPoints1,1) x1(:)];
%             w(countOld:count)=w1(i)*w1(j)*w1(:);
        end
    end  
    
% figure(12)
% plot3(x(:,1),x(:,2),x(:,3),'*')
end