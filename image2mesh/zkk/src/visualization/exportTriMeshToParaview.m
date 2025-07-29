function [] = exportTriMeshToParaview(mesh,options)

    if(nargin==1 || ~isfield(options,'outScreen'))
        outScreen = true;
    else
        outScreen = options.outScreen;
    end
    if(outScreen)
        fprintf('Exporting to paraview...\n')
    end

    exportName = mesh.name;
    
    numNodes = size(mesh.X,1);
    numElems = size(mesh.T,1);
    
    
    nOfFields = 0;
    if(isfield(mesh,'elementField'))
        nOfFields = 1;
    end
    
    if(size(mesh.X,2)==2)
      mesh.X  = [mesh.X zeros(numNodes,1)];
    end   
    
    if(isfield(mesh,'validity'))
        nOfFields = nOfFields+1;
    end
    if(isfield(mesh,'terrain'))
        nOfFields = nOfFields+1;
    end
    
    stream = fopen(['./output/' exportName '.inp'],'w');
    % write header
    fprintf(stream,['%d %d 0 ' int2str(nOfFields) ' 0\n'],numNodes,numElems);
    % write node coordinates
    fprintf(stream,'%d %f %f %f\n', [1:numNodes; mesh.X']);
    % write elements
    fprintf(stream,'%d 0 tri %d %d %d\n', [ 1:numElems; mesh.T']);
    % write visualization properties   
    %printAttributesMesh_paraview(meshLin,stream);
    if(nOfFields==1)
        fprintf(stream,'1 1\n');
        fprintf(stream,'region, None\n');
        fprintf(stream,'%d %f\n',[1:size(mesh.T,1); mesh.elementField']);    
    elseif(nOfFields==2 && isfield(mesh,'validity'))
        fprintf(stream,'2 1 1\n');
        fprintf(stream,'region, None\n');
        fprintf(stream,'validity, None\n');
        fprintf(stream,'%d %d %d\n',[1:size(mesh.T,1); mesh.elementField'; mesh.validity']);   
    elseif(nOfFields==2 && isfield(mesh,'terrain'))
        fprintf(stream,'2 1 1\n');
        fprintf(stream,'region, None\n');
        fprintf(stream,'terrain, None\n');
        fprintf(stream,'%d %d %d\n',[1:size(mesh.T,1); mesh.elementField'; mesh.terrain']);   
    elseif(nOfFields==3)
        fprintf(stream,'3 1 1 1\n');
        fprintf(stream,'region, None\n');
        fprintf(stream,'validity, None\n');
        fprintf(stream,'terrain, None\n');
        fprintf(stream,'%d %d %d %d\n',[1:size(mesh.T,1); mesh.elementField'; mesh.validity'; mesh.terrain']);   
    end
    fclose(stream);

end