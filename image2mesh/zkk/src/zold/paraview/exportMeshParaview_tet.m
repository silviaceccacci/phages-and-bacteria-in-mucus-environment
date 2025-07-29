
function [] = exportMeshParaview_tet(mesh,options)

    % define the CM of each element
    if(isfield(mesh,'CM')==false)
       CM = zeros(size(mesh.X,1),size(mesh.T,1)) ;
       for iElem = 1:size(mesh.T,1)
           CM(:,iElem) = sum(mesh.X(:,mesh.T(iElem,:)),2)/mesh.element.numNod;
       end
       mesh.CM = CM;
    end    
    
    [mesh2D] = getMesh2DFrom3DMesh(mesh);
    exportMeshParaview_tri(mesh2D,options,mesh.element);     
    
    if(isfield(options,'edges') && options.edges)
        exportMeshEdges(mesh,options);
    end
    if(isfield(options,'nodes') && options.nodes)
        exportMeshNodes(mesh,options);
    end
end