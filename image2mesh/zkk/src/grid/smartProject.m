function [M] = smartProject(field, X, T)

    disp('Projecting to mesh (the sophisticated one..)')

    fieldPoints = field.pixelPoints;%[field.x(:) field.y(:)];

    numNodes = size(X,1);
    numElements = size(T,1);
    numFields = size(field.z,3);
    
    elementalField = zeros(numElements,numFields);
    M = zeros(numNodes,numFields);
    countAdjacentElementsToNode = zeros(numNodes,1);
    
    for ielem = 1:numElements
        pixelsInsideElement = findPointsInElement(X(T(ielem,:),:),fieldPoints);
        
        elemPoints = X(T(ielem,:),:);
        numElemPoints = size(elemPoints,1);
        pixelsContainingElemPoints = zeros(numElemPoints,1);
        for ipoin = 1:numElemPoints
            [location] = locate(elemPoints(ipoin,:), field);
            pixelsContainingElemPoints(ipoin) = location.l;
        end
                
%         if(isempty(pixelsInsideElement)) 
%             ielem
%             X(T(ielem,:),:)
%             figure(69); plot( X(T(ielem,[1 2 3 1]),1), X(T(ielem,[1 2 3 1]),2)); grid on;
%             error('Not found field points on element')
%         end
        
        for c=1:numFields
            zi = field.z(:,:,c);
            zi = zi(:);
                            
            elementalField(ielem,c) = ...
                ( sum(zi(pixelsInsideElement)) + sum(zi(pixelsContainingElemPoints)) )/...
                (length(pixelsInsideElement)+length(pixelsContainingElemPoints));
            
            for inode = 1:size(T,2)
                M(T(ielem,inode),c) = M(T(ielem,inode),c) + elementalField(ielem,c);
                if(c==1)
                    countAdjacentElementsToNode(T(ielem,inode)) = countAdjacentElementsToNode(T(ielem,inode))+1;
                end
            end
        end
    end    
    
    %M = bsxfun(@rdivide,M,countAdjacentElementsToNode);
    for inode = 1:numNodes
       if(countAdjacentElementsToNode(inode)==0)
           error('Node info not found');
       end
       M(inode,:) = M(inode,:)/countAdjacentElementsToNode(inode);
    end
    
end