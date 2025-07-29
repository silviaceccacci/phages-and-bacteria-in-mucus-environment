function [isValid,area] = checkMeshValidity(mesh,meshForNormals,minElevFromGround,minArea,lastNonFacadeElem)

    X = mesh.X;
    T = mesh.T;
    Xnormal = meshForNormals.X;
    Tnormal = meshForNormals.T;
    
    numElements = size(T,1);
    isValid = false(numElements,1);
    area = zeros(numElements,1);
    
    for e=1:numElements
       
        if(e>lastNonFacadeElem)% If it is facade
            [isValid(e),area(e)] = checkElemValidity(X(T(e,:),:),Xnormal(Tnormal(e,:),:));
        else
            [isValid(e),area(e)] = checkElemValidity(X(T(e,:),:),[0; 0; 1]);
        end
        
        if(isValid(e))
            if(e>lastNonFacadeElem)% If it is facade
                if(minElevFromGround>0) 
                    isValid(e) = (max(X(T(e,:),3))-min(X(T(e,:),3))) > minElevFromGround;
                end
%                 if(minArea>0)
%                     %isValid(e) = area(e)>minArea;
%                     maxEdge = 0;
%                     elemMinArea = 
%                     elemMinArea = minElevFromGround*minElevFromGround;
%                     isValid(e) = area(e)>minArea;
%                 end        
            end
        end
                
    end

end


function [isValid,area] = checkElemValidity(X,Xnormal)

    if(length(Xnormal(:))==3)
        n = Xnormal;
    else
        v1 = Xnormal(2,:)-Xnormal(1,:);
        v2 = Xnormal(3,:)-Xnormal(1,:);

        e1 = v1'/norm(v1);
        e2 = v2' - (v2*e1)*e1;
        e2 = e2/norm(e2);
        n  = cross(e1,e2);
    end

    v1 = X(2,:)-X(1,:);
    v2 = X(3,:)-X(1,:);
    
    e1 = v1'/norm(v1);
    e2 = v2' - (v2*e1)*e1;
    e2 = e2/norm(e2);
    e3  = cross(e1,e2);
    
    isValid = n'*e3 > 1e-14;
    
    area = norm(cross(v1,v2))/2.0;

%     e1 = Xnormal(1,:)'/norm(Xnormal(1,:));
%     e2 = Xnormal(2,:)' - (Xnormal(2,:)*e1)*e1;
%     e2 = e2/norm(e2);
%     n  = cross(e1,e2);
% 
%     e1 = X(1,:)'/norm(X(1,:));
%     e2 = X(2,:)' - (X(2,:)*e1)*e1;
%     e2 = e2/norm(e2);
%     e3  = cross(e1,e2);
%     
%     isValid = n'*e3 > 1e-14;
      
end
