function [] = exportMeshParaviewSolver(mesh,options)

if isfield(mesh, 'Xp')
    X = mesh.Xp;
    T = mesh.Tp;
else
    X = mesh.X;
    T = mesh.T;
end

if(isfield(options,'name')) 
    exportName = options.name;
    if(isfield(options,'counter'))
        counter    = options.counter;
        exportName = [exportName '_' int2str(counter)];
    end
else
    exportName = 'defaultName';
end

nnodes = size(X,1);
nelems = size(T,1);

nOfFields     = 0;
nOfElemFields = 0;
if(isfield(options,'u'))
    nOfFields = nOfFields +1;
    u=options.u;
else
    nOfFields = nOfFields +1;
    u=zeros(nnodes,2);
end
if(isfield(options,'p'))
    nOfFields = nOfFields +1;
    p=options.p;
else
    nOfFields = nOfFields +1;
    p=zeros(nnodes,1);
end
if(isfield(options,'p0'))
    nOfFields = nOfFields +1;
    p0=options.p0;
else
    nOfFields = nOfFields +1;
    p0=zeros(nnodes,1);
end
if(isfield(options,'k'))
    nOfFields = nOfFields +1;
    k=options.k;
else
    nOfFields = nOfFields +1;
    k=zeros(nnodes,1);
end
if(isfield(options,'m'))
    nOfFields = nOfFields +1;
    m=options.m;
else
    nOfFields = nOfFields +1;
    m=zeros(nnodes,1);
end

sfn = int2str(nOfFields);
sfe = int2str(nOfElemFields);

stream = fopen([exportName '.inp'],'w');
% write header
fprintf(stream,['%d %d ' sfn ' ' sfe ' 0 \n'],nnodes,nelems);
% write node coordinates
fprintf(stream,'%d %e %e 0.0\n', [1:nnodes; X']);
% write elements
fprintf(stream,'%d 0 tri %d %d %d\n', [ 1:nelems; T']);
% write visualization properties   
if nOfFields==2
    fprintf(stream,'2 2 1\n');
    fprintf(stream,'u, None/None\n');
    fprintf(stream,'p, None\n');
    fprintf(stream,'%d %f %f %f \n',[1:nnodes; u'; p']);
elseif nOfFields==3
    fprintf(stream,'3 2 1 1\n');
    fprintf(stream,'u, None/None\n');
    fprintf(stream,'p, None\n');
    fprintf(stream,'k, None\n');
    fprintf(stream,'%d %f %f %f %f\n',[1:nnodes; u'; p'; k']);
elseif nOfFields==4
    fprintf(stream,'4 2 1 1 1\n');
    fprintf(stream,'u, None/None\n');
    fprintf(stream,'p, None\n');
    fprintf(stream,'k, None\n');
    fprintf(stream,'bc, None\n');
    fprintf(stream,'%d %f %f %f %f %f\n',[1:nnodes; u'; p'; k'; m']);
elseif nOfFields==5
    fprintf(stream,'5 2 1 1 1 1\n');
    fprintf(stream,'u, None/None\n');
    fprintf(stream,'p, None\n');
    fprintf(stream,'k, None\n');
    fprintf(stream,'bc, None\n');
    fprintf(stream,'p0, None\n');
    fprintf(stream,'%d %f %f %f %f %f %f\n',[1:nnodes; u'; p'; k'; m'; p0']);
elseif nOfFields==0
    % nothing
else
    error('not implemented nOfFields')
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
    error('not implemented nOfElemFields')
end
    
fclose(stream);
