function [ENglobal]=giveConnectivity_ElementToNode(T,numNodes)

numElements=size(T,1);
numNodesElement=size(T,2);

kTot=numElements*numNodesElement;
i=zeros(kTot,1);
j=zeros(kTot,1);
s=zeros(kTot,1);
k=1;
for iElem=1:numElements
    kn=k+numNodesElement-1;
    i(k:kn)=iElem;
    j(k:kn)=T(iElem,:);
    s(k:kn)=1;
    k=kn+1;
end
ENglobal=sparse(i,j,s,numElements,numNodes);