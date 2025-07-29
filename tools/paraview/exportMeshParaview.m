function [] = exportMeshParaview(X,T,options)
    
    if(nargin<=2) options.nothing = 1; end

	if(isfield(options,'exportName')) 
        exportName = options.exportName;
        if(isfield(options,'counter'))
            counter    = options.counter;
            exportName = [exportName '_' int2str(counter)];
        end
    else
        global iplotg;
        iplotg = iplotg+1;
        exportName = ['defaultName' int2str(iplotg)];
    end
    
    
%     global doSymX
%     if doSymX
%         X(:,1)=-X(:,1);
%     end
%     global doRot90
%     if doRot90
%         xx = X;
%         X(:,1)=xx(:,2);
%         X(:,2)=xx(:,1);
%         clear xx;
%     end
    
    nnodes = size(X,1);
    nelems = size(T,1);
    
    nOfFields     = 0;
    nOfElemFields = 0;
    if(isfield(options,'f'))
        nOfFields = nOfFields +1;
        f=options.f;
    else
        nOfFields = nOfFields +1;
        f=zeros(nnodes,1);
    end
    
    sfn = int2str(nOfFields);
    sfe = int2str(nOfElemFields);
 
    stream = fopen(['./output/paraview/' exportName '.inp'],'w');
    % write header
    fprintf(stream,['%d %d ' sfn ' ' sfe ' 0 \n'],nnodes,nelems);
    % write node coordinates
    fprintf(stream,'%d %f %f 0.0\n', [1:nnodes; X']);
    % write elements
    fprintf(stream,'%d 0 tri %d %d %d\n', [ 1:nelems; T']);
    % write visualization properties   
    if nOfFields==1
        if(size(f,2)==1) 
            fprintf(stream,'1 1\n');
            fprintf(stream,'f, None\n');
            fprintf(stream,'%d %f \n',[1:nnodes; f']);
        else
           fprintf(stream,'1 2\n');
           fprintf(stream,'f, None\n');
           fprintf(stream,'%d %f %f \n',[1:nnodes; f']);
        end

    elseif nOfFields==2
        fprintf(stream,'2 2 1\n');
        fprintf(stream,'u, None/None\n');
        fprintf(stream,'p, None\n');
        fprintf(stream,'%d %f %f %f \n',[1:nnodes; u'; p']);

    elseif nOfFields==6
        fprintf(stream,'6 2 1 1 1 1 1\n');
        fprintf(stream,'u, None/None\n');
        fprintf(stream,'p, None\n');
        fprintf(stream,'f, None\n');
        fprintf(stream,'dist, None\n');
        fprintf(stream,'angle, None\n');
        fprintf(stream,'h, None\n');
        fprintf(stream,'%d %f %f %f %f %f %f %f\n',...
            [1:nnodes; u'; p'; f'; dist'; angle';h']);

    elseif nOfFields==0
        % nothing
    else
        error('not implemented node fields')
    end
    if nOfElemFields==1
        fprintf(stream,'1 1\n');
        fprintf(stream,'d, None\n');
        fprintf(stream,'%d %f\n',[1:nelems; d']);

    elseif nOfElemFields==2
        fprintf(stream,'2 1 2\n');
        fprintf(stream,'d, None\n');
        fprintf(stream,'s, None/None\n');
        fprintf(stream,'%d %f\n',[1:nelems; d'; sourceF']);

    elseif nOfElemFields==0
        % nothing
    else
        error('not implemented elem fields')
    end
        
    fclose(stream);

end