function [mesh] = addMeshConnectivities(mesh)

    [matrixAdjacentElement,matrixLocalFaceAdjacentElement,boundaryNodes,NN]=...
        getMatrixAdjacentElement_general_city(mesh.T,size(mesh.X,1),mesh.element);

    mesh.matrixAdjacentElement = matrixAdjacentElement;
    mesh.matrixLocalFaceAdjacentElement = matrixLocalFaceAdjacentElement;
    %mesh.boundaryNodes = boundaryNodes;
    mesh.boundaryNodesKK = boundaryNodes;
    mesh.NN = NN;
    
    if(length(mesh.boundaryNodesKK) ~= length(mesh.boundaryNodes))
       error('Topological error in 2D mesh: ensure well defined domain') 
    end
end