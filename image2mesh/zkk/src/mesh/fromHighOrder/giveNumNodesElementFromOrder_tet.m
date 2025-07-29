function [ numNod ] = giveNumNodesElementFromOrder_tet(order)

% numNod = 0;
% 
% for p = 0 : order;
% 
%     numNod = numNod + ( (p+1)*(p+2) ) / 2 ;
% 
% end


numNod = round(  ( (order+1)*(order+2)*(order+3) )  /6 );