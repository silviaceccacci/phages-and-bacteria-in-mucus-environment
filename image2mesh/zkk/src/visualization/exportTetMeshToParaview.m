function [] = exportTetMeshToParaview(mesh,options)
    
    fprintf('Exporting to paraview...\n')

    exportName = mesh.name;
    
    numNodes = size(mesh.X,1);
    numElems = size(mesh.T,1);
    
    nOfFields = 0;
    if(isfield(mesh,'quality'))
        nOfFields=nOfFields+1;
    end
        
    stream = fopen(['./output/' exportName '.inp'],'w');
    % write header
    fprintf(stream,['%d %d 0 ' int2str(nOfFields) ' 0\n'],numNodes,numElems);
    % write node coordinates
    fprintf(stream,'%d %f %f %f\n', [1:numNodes; mesh.X']);
    % write elements
    fprintf(stream,'%d 0 tet %d %d %d %d\n', [ 1:numElems; mesh.T']);
    % write visualization properties   
%     %printAttributesMesh_paraview(meshLin,stream);
    if(nOfFields==1)
        fprintf(stream,'1 1\n');
        fprintf(stream,'quality, None\n');
        fprintf(stream,'%d %f\n',[1:size(mesh.T,1); mesh.quality']); 
    end
        
    fclose(stream);

end