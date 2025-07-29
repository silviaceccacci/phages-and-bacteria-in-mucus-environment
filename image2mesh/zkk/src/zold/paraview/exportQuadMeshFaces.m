function [] = exportQuadMeshFaces(mesh,options)
    exportName = options.exportName;
    
    fprintf('Creating the refined quads for the faces \n')
    [meshLin] = defineRefinedLinearMesh(mesh,options);
    nnodesLin = size(meshLin.X,2);
    nelemsLin = size(meshLin.T,1);
 
    if(size(meshLin.X,1)==2)
        meshLin.X = [meshLin.X; zeros(1,size(meshLin.X,2))];
    end
    
    fprintf('Printing the faces \n')
    stream = fopen([exportName '_faces.inp'],'w');
    % write header
    fprintf(stream,'%d %d 0 4 0\n',nnodesLin,nelemsLin);
    % write node coordinates
%     fprintf(stream,'%d %e %e %e\n', [1:nnodesLin; meshLin.X]);
    fprintf(stream,'%d %f %f %f\n', [1:nnodesLin; meshLin.X]);
    % write elements
    fprintf(stream,'%d 0 quad %d %d %d %d \n', [ 1:nelemsLin; meshLin.T']);
    
    
    printAttributesMesh_paraview(meshLin,stream);
    
%     % write visualization properties
%     fprintf(stream,'4 1 1 1 1\n');
%     fprintf(stream,'q, None\n');
%     fprintf(stream,'x, None\n');
%     fprintf(stream,'y, None\n');
%     fprintf(stream,'z, None\n');
%     fprintf(stream,'%d %e %e %e %e \n',[1:nelemsLin; meshLin.quality'; meshLin.CM]);
%     fclose(stream);

end