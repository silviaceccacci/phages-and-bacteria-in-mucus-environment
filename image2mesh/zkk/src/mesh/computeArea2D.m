function [area] = computeArea2D(mesh)

    X = mesh.X;
    T = mesh.T;
    
    numElements = size(T,1);
    area = zeros(numElements,1);
    
    for e=1:numElements
       
        area(e) = computeAreaElem(X(T(e,:),:));

    end

end


function [area] = computeAreaElem(X)

    v1 =[ X(2,:)-X(1,:) 0];
    v2 =[ X(3,:)-X(1,:) 0];
    
    area = norm(cross(v1,v2))/2.0;

end