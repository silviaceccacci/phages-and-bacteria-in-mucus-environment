function [X,T]=removeRepeatedNodes(X,T,numBinsInput)
% by Abel GP
tic;

numNodes = size(X,2);
mapNodes = zeros(numNodes,1);
tolCheck = 1e-12;

if(nargin==2)
    numBins = ceil(sqrt(numNodes));
else
    numBins = numBinsInput;
end


% tic
% for inode=1:numNodes
%     if(mapNodes(inode)==0)
%         listNodes = (inode+1):numNodes;
%         listNodes = listNodes(find(mapNodes(listNodes)==0));
%         vect = bsxfun(@minus,X(:,inode),X(:,listNodes));
%         dist = sqrt(sum(vect.*vect,1));
%         theSameNodes = listNodes(find(dist<tolCheck));
%         mapNodes(theSameNodes) = inode;
%     end
% end
% mapNodes(1:10)
% toc
% 
% mapNodes(:) = 0;

[OT]=OcTree(X','binCapacity',numBins);
for ibin=1:OT.BinCount
    listPointsInBin = find(OT.PointBins==ibin);
    numPointsInBin = length(listPointsInBin);
    for ibinPoint=1:numPointsInBin
        inode = listPointsInBin(ibinPoint);
        if(mapNodes(inode)==0)
            listNodes = listPointsInBin((ibinPoint+1):numPointsInBin);
            listNodes = listNodes(find(mapNodes(listNodes)==0));
            vect = bsxfun(@minus,X(:,inode),X(:,listNodes));
            dist = sqrt(sum(vect.*vect,1));
            theSameNodes = listNodes(find(dist<tolCheck));
            mapNodes(theSameNodes) = inode;
        end
    end
end
% mapNodes(1:10)
% toc
% disp('  ')


nonRepeNodes = find(mapNodes==0);
repeNodes = find(mapNodes);
numNonRepeNodes = length(nonRepeNodes);
newNodeId = zeros(numNodes,1);
newNodeId(nonRepeNodes) = 1:numNonRepeNodes;
newNodeId(repeNodes) = newNodeId(mapNodes(repeNodes));
fprintf('   Number of original nodes: %d\n',numNodes)
fprintf('   Number of non-repeated nodes: %d\n   ',numNonRepeNodes)

X = X(:,nonRepeNodes);
T = newNodeId(T);

toc

% mapNodes(1:10)
% tic
% for inode=1:numNodes
%     if(mapNodes(inode)==0) 
%         for jnode=(inode+1):numNodes
%             if(norm(X(:,inode)-X(:,jnode))<tolCheck)
%                 mapNodes(jnode) = inode;
%             end
%         end
%     end
% end
% toc
% mapNodes(1:10)