function []=computeDiffMeshes(T,X1,X2)

%         L = (max(max(abs(US)))-min(min(abs(US))));
        L=0;
        for iElem = 1:3
            L = max(L,norm(X1(:,T(iElem,1))-X2(:,T(iElem,2))));
            L = max(L,norm(X1(:,T(iElem,2))-X2(:,T(iElem,3))));
            L = max(L,norm(X1(:,T(iElem,3))-X2(:,T(iElem,1))));
        end
        fprintf('mesh %d: abs error %e, rel error %e\n',i,...
            max(norm((X1-X2).^2)),...        
            max(norm((X1-X2).^2))/ L ...
        )

end