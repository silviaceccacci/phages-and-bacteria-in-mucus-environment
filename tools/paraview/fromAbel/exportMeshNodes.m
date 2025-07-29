function [] = exportMeshNodes(mesh,options)
    exportName =options.exportName;
    
    mesh0D.T  = zeros(size(mesh.T,1)*size(mesh.T,2),1);
    mesh0D.X  = zeros(size(mesh.X,1),size(mesh.T,1)*size(mesh.T,2));
    mesh0D.CM = zeros(size(mesh.X,1),size(mesh.T,1)*size(mesh.T,2));
    mesh0D.quality  = zeros(size(mesh0D.T,1),1);
    if(isfield(mesh,'qualities'))
       namesQ = fieldnames(mesh.qualities);
       for i=1:length(namesQ)
           qualitiesSubelements.(namesQ{i}) = zeros(size(mesh0D.quality));
       end
    end
    if(isfield(mesh,'distortions'))
       namesQ = fieldnames(mesh.distortions);
       for i=1:length(namesQ)
           distortionsSubelements.(namesQ{i}) = zeros(size(mesh0D.quality));
       end
    end
    if(isfield(mesh,'errorDisc'))
       namesE = fieldnames(mesh.errorDisc);
       for i=1:length(namesE)
           errorDisc.(namesE{i}) = zeros(size(mesh0D.quality));
       end
    end
    countNodes = 0;
    for iElem = 1:size(mesh.T,1)
        for iNode = 1:size(mesh.T,2)
            countNodes = countNodes +1;
            mesh0D.T(countNodes,1)  = countNodes;
            mesh0D.X(:,countNodes)  = mesh.X(:,mesh.T(iElem,iNode));
            mesh0D.CM(:,countNodes) = mesh.CM(:,iElem);
            if(isfield(mesh,'quality'))
                mesh0D.quality(countNodes) = mesh.quality(iElem);
            else
                mesh0D.quality(countNodes) = 1;
            end
           if(isfield(mesh,'qualities'))
               for i=1:length(namesQ)
                   qualitiesSubelements.(namesQ{i})(countNodes) = mesh.qualities.(namesQ{i})(iElem);
               end
           end
           if(isfield(mesh,'distortions'))
               for i=1:length(namesQ)
                   distortionsSubelements.(namesQ{i})(countNodes) = mesh.distortions.(namesQ{i})(iElem);
               end
           end
           if(isfield(mesh,'errorDisc'))
               for i=1:length(namesE)
                   errorDisc.(namesE{i})(countNodes) = mesh.errorDisc.(namesE{i})(iElem);
               end
           end
        end
    end
    if(isfield(mesh,'qualities'))
       for i=1:length(namesQ)
           mesh0D.qualities.(namesQ{i}) = qualitiesSubelements.(namesQ{i})(1:countNodes);
       end
    end
    if(isfield(mesh,'distortions'))
       for i=1:length(namesQ)
           mesh0D.distortions.(namesQ{i}) = distortionsSubelements.(namesQ{i})(1:countNodes);
       end
    end
    if(isfield(mesh,'errorDisc'))
       for i=1:length(namesE)
           mesh0D.errorDisc.(namesE{i}) = errorDisc.(namesE{i})(1:countNodes);
       end
    end
    
    fprintf('Printing the nodes \n')
    stream = fopen([exportName '_nodes.inp'],'w');
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
    fprintf(stream,['%d %d 0 ' int2str(nOfFields) ' 0\n'],size(mesh0D.X,2),size(mesh0D.T,1));
%     fprintf(stream,'%d %d 0 3 0\n',size(mesh0D.X,2),size(mesh0D.T,1));
    % write node coordinates
    fprintf(stream,'%d %f %f %f\n', [1:size(mesh0D.X,2); mesh0D.X]);
    % write elements
    fprintf(stream,'%d 0 pt %d\n', [ 1:size(mesh0D.T,1); mesh0D.T']);
    % write visualization properties
    printAttributesMesh_paraview(mesh0D,stream);
%     fprintf(stream,'3 1 1 1\n');
%     fprintf(stream,'x, None\n');
%     fprintf(stream,'y, None\n');
%     fprintf(stream,'z, None\n');
%     fprintf(stream,'%d %e %e %e \n',[1:size(mesh0D.T,1); mesh0D.CM]);
    fclose(stream);

end


%     if(isfield(mesh,'qualities'))
%        qualitiesSubelements.IdGeo =  zeros(size(mesh0D.quality));
%        qualitiesSubelements.IdIni =  zeros(size(mesh0D.quality));
%        qualitiesSubelements.ScJac =  zeros(size(mesh0D.quality));
%     end
%             if(isfield(mesh,'qualities'))
%                qualitiesSubelements.IdGeo(countNodes) =  mesh.qualities.IdGeo(iElem);
%                qualitiesSubelements.IdIni(countNodes) =  mesh.qualities.IdIni(iElem);
%                qualitiesSubelements.ScJac(countNodes) =  mesh.qualities.ScJac(iElem);
%             end