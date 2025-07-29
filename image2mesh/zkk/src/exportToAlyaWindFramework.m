function []=exportToAlyaWindFramework(fileName,exportToAlya,...
    mesh3D,mesh_facade3D,mesh2D,...
    fields,field_topo,translationVector_saved,domainLimitsIni)

if(exportToAlya)
    tic; fprintf('Writting Alya files...\n')
    mesh3D.fields = fields;
    mesh3D.domainLimits = domainLimitsIni;

    fileName = [fileName '.WindMesh'];

    exportMetaData(fileName,mesh2D,translationVector_saved)

    [meshTopo] = exportCuts(fileName,field_topo,mesh2D,mesh_facade3D,false,translationVector_saved);

    options.meshTopo = meshTopo;
    options.meshFacades = mesh_facade3D;
    exportAlya(mesh3D,fileName,options);
end

% mesh3D = rmfield(mesh3D,'fields');
% mesh3D = rmfield(mesh3D,'domainLimits');
% mesh3D = rmfield(mesh3D,'boundaryNodes');
% if(isfield(mesh3D,'nodesSmooth'))
%     mesh3D = rmfield(mesh3D,'nodesSmooth');
% end
% if(isfield(mesh3D,'idealElements'))
%     mesh3D = rmfield(mesh3D,'idealElements');
% end

