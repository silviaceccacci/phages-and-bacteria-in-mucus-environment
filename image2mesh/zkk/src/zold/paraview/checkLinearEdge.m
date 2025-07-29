function [logic] = checkLinearEdge(X,options)

    logic=true;
    
    if(nargin>1)
        tol = options.tol;
    else
        tol = 1e-6; %1e-4;%1e-6;
    end
    
    numNod= size(X,2);
    e = X(:,numNod)-X(:,1);
    for iPiece = 2:(numNod-1)
        ei = X(:,iPiece)-X(:,1);
        if(size(X,1)==2)
            logic = find(abs(det([e ei]))<tol);
        elseif(size(X,1)==3)
%             eaux1 = [1;0;0];
%             eaux2 = [0;1;0];
%             eaux3 = [0;0;1];
%             logic = min([ (abs(det([e ei eaux1]))<tol) ...
%                           (abs(det([e ei eaux2]))<tol) ...
%                           (abs(det([e ei eaux3]))<tol)]);
            logic = min( [ ...
                    (abs(det([e([1 2]) ei([1 2])]))<tol)...
                    (abs(det([e([2 3]) ei([2 3])]))<tol)...
                    (abs(det([e([1 3]) ei([1 3])]))<tol) ] );
        end
        if(logic==false)
            break;
        end
    end
    
end