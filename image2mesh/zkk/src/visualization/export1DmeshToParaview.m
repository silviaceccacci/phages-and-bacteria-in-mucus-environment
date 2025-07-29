function [] = export1DmeshToParaview(mesh,exportName,options)

    if(nargin==2 || ~isfield(options,'outScreen'))
        outScreen = true;
    else
        outScreen = options.outScreen;
    end
    if(outScreen)
        fprintf('Exporting to paraview...\n')
    end
    
    numNodes = size(mesh.X,1);
    numElems = size(mesh.T,1);
    
    nOfFields = 0;%1;
    
    if(size(mesh.X,2)==2)
      mesh.X  = [mesh.X zeros(numNodes,1)];
    end   
    
%     if(isfield(mesh,'validity'))
%         nOfFields = nOfFields+1;
%     end
    
    stream = fopen(['./output/' exportName '.inp'],'w');
    % write header
    fprintf(stream,['%d %d 0 ' int2str(nOfFields) ' 0\n'],numNodes,numElems);
    % write node coordinates
    fprintf(stream,'%d %f %f %f\n', [1:numNodes; mesh.X']);
    % write elements
    fprintf(stream,'%d 0 line %d %d\n', [ 1:numElems; mesh.T']);
    % write visualization properties   
    %printAttributesMesh_paraview(meshLin,stream);
%     if(nOfFields==1)
%         fprintf(stream,'1 1\n');
%         fprintf(stream,'region, None\n');
%         fprintf(stream,'%d %f\n',[1:size(mesh.T,1); mesh.elementField']);    
%     elseif(nOfFields==2)
%         fprintf(stream,'2 1 1\n');
%         fprintf(stream,'region, None\n');
%         fprintf(stream,'validity, None\n');
%         fprintf(stream,'%d %f %d\n',[1:size(mesh.T,1); mesh.elementField'; mesh.validity']);   
%     end
    fclose(stream);

end