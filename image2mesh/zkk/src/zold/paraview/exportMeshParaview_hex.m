function [] = exportMeshParaview_hex(mesh,options)
    [mesh2D] = getMesh2DFrom3DMesh(mesh);
    exportMeshParaview_quad(mesh2D,options);    
end