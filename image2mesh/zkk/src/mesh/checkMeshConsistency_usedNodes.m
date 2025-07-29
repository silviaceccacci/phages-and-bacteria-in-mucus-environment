function [consistency,mesh]=checkMeshConsistency_usedNodes(mesh)

    consistency = size(mesh.X,1)==length(unique(mesh.T));

    if(~consistency)
        usedNodes = unique(mesh.T);
        mapToNew = zeros(size(mesh.X,1),1);
        mapToNew(usedNodes) = 1:length(usedNodes);
        
        mesh.X = mesh.X(usedNodes,:);
        mesh.T = mapToNew(mesh.T); 
    end

%     if(~consistency)
%         disp('   Nodes not present in topology')
%     else
%         disp('   Consistent mesh')
%     end
    
end