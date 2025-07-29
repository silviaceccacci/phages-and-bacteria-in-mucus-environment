function []=exportBadQualityMesh(mesh0,mesh1,meshp1)

    badQuality = 0.1;

    deg = find(mesh0.quality<badQuality);

    mesh0.idealGeometric = mesh1.idealGeometric;
    mesh0.idealInitial   = mesh1.idealInitial;
    mesh0.idealElements  = mesh1.idealElements;

    exportSubmesh(mesh0,deg);
    exportSubmesh(mesh1,deg);
    
    if(nargin>2)
        meshp1.idealGeometric = mesh1.idealGeometric;
        meshp1.idealInitial   = mesh1.idealInitial;
        meshp1.idealElements  = mesh1.idealElements;
        exportSubmesh(meshp1,deg);
    end
end

function []=exportSubmesh(mesh,deg)


    meshDeg = mesh;
    meshDeg.T = mesh.T(deg,:); 
    meshDeg.quality = mesh.quality(deg);
    meshDeg.idealGeometric = mesh.idealGeometric(deg,:,:);
    meshDeg.idealInitial   = mesh.idealInitial(deg,:,:);
    meshDeg.idealElements  = mesh.idealElements(deg,:,:);
    if( isfield(mesh,'qualities') ) 
        meshDeg.qualities.IdGeo = mesh.qualities.IdGeo(deg);
        meshDeg.qualities.IdIni = mesh.qualities.IdIni(deg);
        meshDeg.qualities.ScJac = mesh.qualities.ScJac(deg);
    end
    
    if( isfield(mesh,'distortions') )
        meshDeg.distortions.IdGeo = mesh.distortions.IdGeo(deg);
        meshDeg.distortions.IdIni = mesh.distortions.IdIni(deg);
        meshDeg.distortions.ScJac = mesh.distortions.ScJac(deg);
    end
    
    meshDeg.fileName = [mesh.fileName 'Deg'];
    
    

    options.exportOrder = 4;
    options.edges = 1;
    options.nodes = 0;
    options.solid = 1;
    options.reduction = 0;
    
    exportMeshParaview(meshDeg,options);

end


