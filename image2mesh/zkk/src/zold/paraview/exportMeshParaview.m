function [] = exportMeshParaview(mesh,options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exports a .inp file of the mesh
% Options:
%    -options.reduction = 1 -> export as linear the straight-sided faces
%    -options.exportOrder   -> exports the mesh with subelements as if were of exportOrder degree
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(isfield(mesh,'distortions'))
    mesh = rmfield(mesh,'distortions');
end

    rute = 'output/';
    if(isfield(mesh,'fileName'))
        exportName = [rute mesh.fileName];
    else
        exportName = [rute 'unknownName'];
    end
    options.exportName = exportName;

    if(isfield(options,'reduction'))
        options.activeCostReduction = options.reduction;
    else
        options.activeCostReduction = true;
    end

    if(isfield(options,'exportOrder')==false)
        exportOrder = 6;
        options.exportOrder = max(exportOrder,mesh.element.order);
    end
    
    fprintf('Exporting inp file: %s\n',options.exportName);

    switch mesh.element.type
        case 'tri'
            exportMeshParaview_tri(mesh,options);
        case 'tet'
            if(nargin>1 && isfield(options,'solid') && options.solid)       
                exportMeshParaview_tetSolid(mesh,options);
            else
                exportMeshParaview_tet(mesh,options);
            end
        case 'quad'
            exportMeshParaview_quad(mesh,options);
        case 'hex'
            exportMeshParaview_hex(mesh,options);
        otherwise
            fprintf('element not implemented \n')
    end

end


% use  
%       names = fieldnames(aStruct)
% to extract the name of the fields of a struct 
% and then be able to print them



















