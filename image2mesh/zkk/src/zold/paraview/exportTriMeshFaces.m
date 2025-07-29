function [] = exportTriMeshFaces(mesh,options)
    exportName = options.exportName;
    
    fprintf('Creating the refined triangulation for the faces \n')
    [meshLin] = defineRefinedLinearMesh(mesh,options);
    nnodesLin = size(meshLin.X,2);
    nelemsLin = size(meshLin.T,1);
 
    fprintf('Printing the faces \n')
    stream = fopen([exportName '_faces.inp'],'w');
    % write header
    nOfFields = 3;
    if(isfield(mesh,'quality'))
        nOfFields = nOfFields +1;
    end
    if(isfield(mesh,'qualities'))
        nOfFields = nOfFields +3;
    end
    if(isfield(mesh,'distortions') ||  isfield(mesh,'errorDisc'))
        nOfFields = nOfFields +3;
    end
    fprintf(stream,['%d %d 0 ' int2str(nOfFields) ' 0\n'],nnodesLin,nelemsLin);
    % write node coordinates
%     fprintf(stream,'%d %e %e %e\n', [1:nnodesLin; meshLin.X]);
    fprintf(stream,'%d %f %f %f\n', [1:nnodesLin; meshLin.X]);
    % write elements
    fprintf(stream,'%d 0 tri %d %d %d\n', [ 1:nelemsLin; meshLin.T']);
    % write visualization properties   
    printAttributesMesh_paraview(meshLin,stream);
%     if(isfield(meshLin,'qualities')==false)
%         fprintf(stream,'4 1 1 1 1\n');
%         fprintf(stream,'q, None\n');
%         fprintf(stream,'x, None\n');
%         fprintf(stream,'y, None\n');
%         fprintf(stream,'z, None\n');
%         fprintf(stream,'%d %e %e %e %e \n',[1:nelemsLin; meshLin.quality'; meshLin.CM]);
%     else
%         fprintf(stream,'7 1 1 1 1 1 1 1\n');
%         fprintf(stream,'q, None\n');
%         fprintf(stream,'x, None\n');
%         fprintf(stream,'y, None\n');
%         fprintf(stream,'z, None\n');
%         fprintf(stream,'IdGeo, None\n');
%         fprintf(stream,'IdIni, None\n');
%         fprintf(stream,'ScJac, None\n');
%         fprintf(stream,'%d %e %e %e %e %e %e %e \n',...
%             [1:nelemsLin; meshLin.quality'; meshLin.CM;...
%             meshLin.qualities.IdGeo';...
%             meshLin.qualities.IdIni';...
%             meshLin.qualities.ScJac']);
%     end     
        
        
        
    fclose(stream);

end