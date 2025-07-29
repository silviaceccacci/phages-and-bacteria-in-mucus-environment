function [] = exportMeshEdges(mesh,options)
    fprintf('Printing the edges \n')
    exportName = options.exportName;
%     fprintf('Creating the refined lines for the edges \n')
    [mesh1D] = getEdgeCurves(mesh,options);
    nnodes1D = size(mesh1D.X,2);
    nelems1D = size(mesh1D.T,1);
    
    stream = fopen([exportName '_edges.inp'],'w');   
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
    fprintf(stream,['%d %d 0 ' int2str(nOfFields) ' 0\n'],nnodes1D,nelems1D);
%     fprintf(stream,'%d %d 0 4 0\n',nnodes1D,nelems1D);
    % write node coordinates
    fprintf(stream,'%d %f %f %f\n', [1:nnodes1D; mesh1D.X]);
    % write elements
    fprintf(stream,'%d 0 line %d %d\n', [ 1:nelems1D; mesh1D.T']);
    % write visualization properties
    printAttributesMesh_paraview(mesh1D,stream);
%     fprintf(stream,'4 1 1 1 1\n');
%     fprintf(stream,'q, None\n');
%     fprintf(stream,'x, None\n');
%     fprintf(stream,'y, None\n');
%     fprintf(stream,'z, None\n');
%     fprintf(stream,'%d %e %e %e %e \n',[1:nelems1D; mesh1D.quality'; mesh1D.CM]);
    fclose(stream);
end