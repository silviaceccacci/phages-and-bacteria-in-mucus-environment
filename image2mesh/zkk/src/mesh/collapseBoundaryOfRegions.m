function [polylines,doCollapse,doCloseRegions,countLoops] =...
            collapseBoundaryOfRegions(mesh,tolCollapse,tolClose,doCollapse,doCloseRegions,...
            countLoops,maxLoops,doCheckArea,minElementAreaAllowed)

isCollapsed = false;
isClosed = false;

polylines.X = [];
polylines.T = [];

numEliminatedRegions = 0;
numCollapsedPoints = 0;
numClosedRegions = 0;

numRegions = size(mesh.fieldElements,1);
for iregion = (mesh.groundRegion+1):numRegions
    regionElements = mesh.fieldElements{iregion};
    
    regionElementsArea = mesh.area(regionElements);
    
    maxRegionElementsArea = max(regionElementsArea);
    minRegionElementsArea = min(regionElementsArea);
    if(doCheckArea && minRegionElementsArea<minElementAreaAllowed)
            numEliminatedRegions = numEliminatedRegions+1;
    elseif(maxRegionElementsArea<=tolCollapse^2)
        numEliminatedRegions = numEliminatedRegions+1;
    else
        %% construct boundary loop
        adjElem = mesh.matrixAdjacentElement(regionElements,:);
        linEdgeIds = find(mesh.elementField(adjElem)~=iregion);
        [bel,bed] = ind2sub(size(adjElem),linEdgeIds);

        edgeNodes = [1 2; 2 3; 3 1];
        bEdgeNodes = zeros(length(bel),2);
        for iaux = 1:length(bel)
            bEdgeNodes((iaux),1) = mesh.T(regionElements(bel(iaux)),edgeNodes(bed(iaux),1));    
            bEdgeNodes((iaux),2) = mesh.T(regionElements(bel(iaux)),edgeNodes(bed(iaux),2)); 
        end
        
        numBEdges = size(bEdgeNodes,1);
        remainingE = 2:numBEdges;
        bnodes = zeros(numBEdges,1);%zeros(numBEdges+1,1);
        bnodes(1:2) = bEdgeNodes(1,:) ;

        countNod = 2;
        while(~isempty(remainingE))
            selectedEdge = (find(bnodes(countNod)==bEdgeNodes(remainingE,1)));

            nextEdge = remainingE(selectedEdge);

            countNod = countNod+1;
            bnodes(countNod) = bEdgeNodes(nextEdge,2);

            remainingE(selectedEdge) = [];
        end
        if(doCollapse)
            %% collapse boundary loop
            newb = [];
            count = 0;
            while(length(newb)~=length(bnodes))
                if(count > 0)
                    bnodes = newb;
                end
                newb = zeros(size(bnodes));
                countNewNodes = 0;
                prevCheck = 1;
                for inode=2:length(newb)
                    if( norm(mesh.X(bnodes(inode),:)-mesh.X(bnodes(prevCheck),:)) > tolCollapse )
                        countNewNodes=countNewNodes+1;
                        newb(countNewNodes) = bnodes(inode);
                        prevCheck = inode;
                    else
                        %mesh.X(bnodes(inode+1),:) = ( mesh.X(bnodes(inode),:)+mesh.X(bnodes(inode+1),:) )/2.0;
                        numCollapsedPoints = numCollapsedPoints+1;
                        isCollapsed = true;
                    end
                end
                newb = [newb(countNewNodes); newb(1:countNewNodes)];
                count = count+1;
            end
            %% configure new polyline to mesh
            newb = newb(1:(end-1));
            if(length(newb)>1)
                lastNewNode = size(polylines.X,1);
                newIndexNodes = [ (lastNewNode+1):(lastNewNode+length(newb)) ]';
                polylines.T = [polylines.T
                    newIndexNodes(1:(end-1)) newIndexNodes(2:end) 
                    newIndexNodes(end)       newIndexNodes(1)   ];
%                 polylines.T = [polylines.T
%                     ((size(polylines.X,1)+1):(size(polylines.X,1)+length(newb)-1))' ,...
%                     ((size(polylines.X,1)+2):(size(polylines.X,1)+length(newb))  )' 
%                     (size(polylines.X,1)+length(newb))  (size(polylines.X,1)+1)   ];
                polylines.X = [polylines.X; mesh.X(newb,:)];
            end
        elseif(doCloseRegions)
            bnodes = bnodes(1:(end-1));
            newPair = [];
            for inode = 1:length(bnodes)
                % find nodes close under the tolerance
                searchNodes = (inode+2):length(bnodes);
                searchNodes = searchNodes(find(mesh.NN(bnodes(searchNodes),bnodes(inode))));
                diff = bsxfun(@minus,mesh.X(bnodes(searchNodes),:),mesh.X(bnodes(inode),:));
                dist = sqrt(sum(diff.^2,2));
                [mind,ind] = min(dist);
                if(mind<tolClose)
                    closestNode = searchNodes(ind);
%                     if(mesh.NN(bnodes(inode),bnodes(closestNode)))
                        %newPair = [newPair; bnodes(inode) bnodes(closestNode)];
                        newPair = [newPair; inode closestNode];
                        %break;
%                     else
%                         dist(ind) = inf;
%                         [mind,ind] = min(dist);
%                     end
                end
            end
            % add to the polylines the contour polylines and the most small
            % non boundary edge under the tolerance (or maybe the first one
            % to make it faster
            lastNewNode = size(polylines.X,1);
            newIndexNodes = [ (lastNewNode+1):(lastNewNode+length(bnodes)) ]';
            polylines.T = [polylines.T
                newIndexNodes(1:(end-1)) newIndexNodes(2:end) 
                newIndexNodes(end)       newIndexNodes(1)   ];
            polylines.X = [polylines.X; mesh.X(bnodes,:)];
            if(~isempty(newPair))
                polylines.T = [polylines.T
                    lastNewNode+newPair   ];
                numClosedRegions = numClosedRegions+1;
                isClosed = true;
            end
        elseif(doCheckArea)
            %Add everything
            lastNewNode = size(polylines.X,1);
            newIndexNodes = [ (lastNewNode+1):(lastNewNode+length(bnodes)) ]';
            polylines.T = [polylines.T
                newIndexNodes(1:(end-1)) newIndexNodes(2:end) 
                newIndexNodes(end)       newIndexNodes(1)   ];
            polylines.X = [polylines.X; mesh.X(bnodes,:)];
        else
            error('non valid option');
        end
    end
end

if(doCollapse)
    fprintf('      Num removed regions: %d\n',numEliminatedRegions)
elseif(doCloseRegions)
    fprintf('      Num closed regions: %d\n',numClosedRegions) 
elseif(doCheckArea)
    fprintf('      Num removed regions: %d\n',numEliminatedRegions)
end


% if(doCollapse)
%     if(doCloseRegions)
%         doCollapse = isCollapsed;
%     else
%         doCollapse = isCollapsed; 
%         doCloseRegions = ~isCollapsed;
%     end
% %     doCollapse = isCollapsed;
% %     if(doCloseRegions && ~doCollapse)
% %         doCloseRegions = false;
% %     else
% %         doCloseRegions = ~isCollapsed;
% %     end
% %     doCloseRegions = ~isCollapsed;
% elseif(doCloseRegions)
%     doCollapse = isClosed;
%     doCloseRegions = true;%false;
% end

if(doCollapse && doCloseRegions)
    countLoops = countLoops+1;
    if(countLoops>=maxLoops)
        doCollapse = false;
        doCloseRegions = false;
    else
        doCollapse = isCollapsed;
        doCloseRegions = false;
    end
else
    if(doCollapse)
        doCollapse = isCollapsed;
        doCloseRegions = ~isCollapsed;
    elseif(doCloseRegions)
        if(isClosed)
            doCollapse = true;
            doCloseRegions = true;
        else
            doCollapse = false;
            doCloseRegions = false;
        end
    end
end


end
