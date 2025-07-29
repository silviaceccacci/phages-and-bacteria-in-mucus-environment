function []=printAttributesMesh_paraview(mesh,stream)

    % write visualization properties   
    if(isfield(mesh,'qualities')==true && isfield(mesh,'errorDisc')==true)
        fprintf(stream,'10 1 1 1 1 1 1 1 1 1 1\n');
        fprintf(stream,'q, None\n');
        fprintf(stream,'x, None\n');
        fprintf(stream,'y, None\n');
        fprintf(stream,'z, None\n');
        fprintf(stream,'QIdGeo, None\n');
        fprintf(stream,'QIdIni, None\n');
        fprintf(stream,'QScJac, None\n');
        fprintf(stream,'EDiscRef, None\n');
        fprintf(stream,'EDiscPhys, None\n');
        fprintf(stream,'EDiscMax, None\n');
        fprintf(stream,'%d %f %f %f %f %f %f %f %f %f %f\n',...
            [1:size(mesh.T,1); mesh.quality'; mesh.CM;...
            mesh.qualities.IdGeo';...
            mesh.qualities.IdIni';...
            mesh.qualities.ScJac';...
            mesh.errorDisc.eRef';...
            mesh.errorDisc.ePhy';...
            mesh.errorDisc.eMax']);
    elseif( isfield(mesh,'qualities')==true &&  isfield(mesh,'distortions')==true )
        fprintf(stream,'10 1 1 1 1 1 1 1 1 1 1\n');
        fprintf(stream,'q, None\n');
        fprintf(stream,'x, None\n');
        fprintf(stream,'y, None\n');
        fprintf(stream,'z, None\n');
        fprintf(stream,'QIdGeo, None\n');
        fprintf(stream,'QIdIni, None\n');
        fprintf(stream,'QScJac, None\n');
        fprintf(stream,'DIdGeo, None\n');
        fprintf(stream,'DIdIni, None\n');
        fprintf(stream,'DScJac, None\n');
        fprintf(stream,'%d %f %f %f %f %f %f %f %f %f %f\n',...
            [1:size(mesh.T,1); mesh.quality'; mesh.CM;...
            mesh.qualities.IdGeo';...
            mesh.qualities.IdIni';...
            mesh.qualities.ScJac';
            mesh.distortions.IdGeo';...
            mesh.distortions.IdIni';...
            mesh.distortions.ScJac']);
    elseif( isfield(mesh,'qualities')==true )
        fprintf(stream,'7 1 1 1 1 1 1 1\n');
        fprintf(stream,'q, None\n');
        fprintf(stream,'x, None\n');
        fprintf(stream,'y, None\n');
        fprintf(stream,'z, None\n');
        fprintf(stream,'QIdGeo, None\n');
        fprintf(stream,'QIdIni, None\n');
        fprintf(stream,'QScJac, None\n');
        fprintf(stream,'%d %f %f %f %f %f %f %f\n',...
            [1:size(mesh.T,1); mesh.quality'; mesh.CM;...
            mesh.qualities.IdGeo';...
            mesh.qualities.IdIni';...
            mesh.qualities.ScJac']);
    elseif( isfield(mesh,'errorDisc')==true )
        fprintf(stream,'7 1 1 1 1 1 1 1\n');
        fprintf(stream,'q, None\n');
        fprintf(stream,'x, None\n');
        fprintf(stream,'y, None\n');
        fprintf(stream,'z, None\n');
        fprintf(stream,'EDiscRef, None\n');
        fprintf(stream,'EDiscPhys, None\n');
        fprintf(stream,'EDiscMax, None\n');
        fprintf(stream,'%d %f %f %f %f %f %f %f\n',...
            [1:size(mesh.T,1); mesh.quality'; mesh.CM;...
            mesh.errorDisc.eRef';...
            mesh.errorDisc.ePhy';...
            mesh.errorDisc.eMax']);
    elseif(isfield(mesh,'quality')==true && isfield(mesh,'qualities')==false && isfield(mesh,'errorDisc')==false)
        fprintf(stream,'4 1 1 1 1\n');
        fprintf(stream,'q, None\n');
        fprintf(stream,'x, None\n');
        fprintf(stream,'y, None\n');
        fprintf(stream,'z, None\n');
        fprintf(stream,'%d %f %f %f %f\n',[1:size(mesh.T,1); mesh.quality'; mesh.CM;]);
    elseif(isfield(mesh,'quality')==false && isfield(mesh,'qualities')==false && isfield(mesh,'errorDisc')==false && isfield(mesh,'CM')==true)
        fprintf(stream,'3 1 1 1\n');
        fprintf(stream,'x, None\n');
        fprintf(stream,'y, None\n');
        fprintf(stream,'z, None\n');
        fprintf(stream,'%d %f %f %f\n',[1:size(mesh.T,1); mesh.CM;]);
    end     
    
end