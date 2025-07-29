function [insidePoins] = findPointsInElement(Xelem,points)

    nextNode = [2 3 1];
    edges = [ 1 2 ; 2 3 ; 3 1];
    isInside = ones(size(points,1),1);
    
    numEdges = size(edges,1);
    for iedge=1:numEdges
       
        n1 = edges(iedge,1);
        n2 = edges(iedge,2);
        n3 = nextNode(n2);
        vedge = Xelem(n2,:) - Xelem(n1,:);
        normals = [ -vedge(2),  vedge(1) 
                     vedge(2), -vedge(1) ];

        normalId =  find(normals*(Xelem(n3,:) - Xelem(n1,:))'>0); 
        normal = normals(normalId,:)';
        
        n = bsxfun(@minus,points,Xelem(n1,:));

        %isInside = isInside && (sign(n*normal)==1);%( (vedge*normal) > -1e-14);  
        isInside = min(isInside,sign(n*normal)>-1);
    end

    insidePoins = find(isInside);
end
