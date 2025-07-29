function [] = exportMeshParaview_tetSolid(mesh,options)

     exportName = options.exportName;
    if(isfield(options,'plotOrder'))
        options.exportOrder = options.plotOrder;
    end
    options
    
    if(isfield(mesh,'CM')==false)
       CM = zeros(size(mesh.X,1),size(mesh.T,1)) ;
       for iElem = 1:size(mesh.T,1)
           CM(:,iElem) = sum(mesh.X(:,mesh.T(iElem,:)),2)/mesh.element.numNod;
       end
       mesh.CM = CM;
    end    

    [nodalDistortion,meshLin] = computeNodalDisortion(options,mesh);
    
    nnodesLin = size(meshLin.X,2);
    nelemsLin = size(meshLin.T,1);
 
    fprintf('Printing the solid tets\n')
    stream = fopen([exportName '_tetSolid.inp'],'w');
    
    % write header
    writeHeader(stream,meshLin,nnodesLin,nelemsLin);
    % write node coordinates
    fprintf(stream,'%d %f %f %f\n', [1:nnodesLin; meshLin.X]);
    % write elements
    switch meshLin.element.type
        case 'tet'
        fprintf(stream,'%d 0 tet %d %d %d %d\n', [ 1:nelemsLin; meshLin.T']);
        case 'tri'
        fprintf(stream,'%d 0 tri %d %d %d\n', [ 1:nelemsLin; meshLin.T']);
            
    end
    indexos = find(isinf(nodalDistortion));
    nodalDistortion(indexos) = 1e6;
    printNodeAttributes_paraview(stream,meshLin,nodalDistortion);
    printAttributesMesh_paraview(meshLin,stream);  

    fclose(stream);
    
    %%
    if(isfield(options,'edges') && options.edges)
        exportMeshEdges(mesh,options);
    end
    if(isfield(options,'nodes') && options.nodes)
        exportMeshNodes(mesh,options);
    end
end

function []=writeHeader(stream,mesh,nnodesLin,nelemsLin)

    nOfNodeFields = 2;

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
    fprintf(stream,['%d %d ' int2str(nOfNodeFields) ' ' int2str(nOfFields) ' 0\n'],...
        nnodesLin,nelemsLin);

end














