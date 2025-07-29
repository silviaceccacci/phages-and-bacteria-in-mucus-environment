function [] = exportSurfaceInterest(mesh,interest)

    [meshSurface]=defineSurfaceInterestMesh(mesh,interest);

    meshSurface.fileName = [mesh.fileName '_object'];

    options.reduction = false;
    options.factor = 3;
    options.discontinuous = 0;
    options.edges=1;
    options.exportOrder = 6;
    exportMeshParaview(meshSurface,options);

end

 