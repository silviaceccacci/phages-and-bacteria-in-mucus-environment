function [] = exportMeshParaview_quad(mesh,options)

    nnodes = size(mesh.X,2);
    nelems = size(mesh.T,1);  
    element = mesh.element;   
    
    % define the CM of each element
    if(isfield(mesh,'CM')==false)
       CM = zeros(size(mesh.X,1),nelems) ;
       for iElem = 1:nelems
           CM(:,iElem) = sum(mesh.X(:,mesh.T(iElem,:)),2)/element.numNod;
       end
       mesh.CM = CM;
    end
        
    %% write faces and edges
    exportQuadMeshFaces(mesh,options);
    if(isfield(options,'edges') && options.edges)
        exportMeshEdges(mesh,options);
    end
    if(isfield(options,'nodes') && options.nodes)
        exportMeshNodes(mesh,options);
    end
        
end