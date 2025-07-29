function [polylines,fields,translationVector,translationVector_saved]...
    =readInputCityMesh(...
    inpName,geomType,fileName_point,fileName,...
    lidarFieldName,topoFieldName,field_lidar,field_topo,...
    exportStepsToParaview,exportLidarParaview,exportTopoParaview)

tic; fprintf('Reading...\n')
[polylines,translationVector]=read1DGeometry(inpName,geomType,fileName_point);
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

% polylines.X = polylines.X/1000;
% translationVector = translationVector/1000;
% warning('----> XXXXXX canviant a metres de manera manualllllllllll XXXXX <<<<--------------------------')


%translationVector(:) = 0;
%warning('He canviat el translation point a zero!!!! ARA VA?????')

if(exportStepsToParaview)
	export1DmeshToParaview(polylines,[fileName '_cadastre'])
    %error('para')
end

% % to get just some of the buildings:
%polylines.T = polylines.T(1:200,:); newNodes = unique(polylines.T(:)); polylines.X = polylines.X(newNodes,:);
%mapOldToNew = zeros(size(polylines.X,1),1); mapOldToNew(newNodes) = 1:length(newNodes); polylines.T(:,:) = mapOldToNew(polylines.T);
    
if(~isempty(lidarFieldName))
    disp('Read lidar field')
%     load BCN_lidarClipped_field;
%     %load(lidarFieldName)     
%     field_lidar = field;
    if(exportLidarParaview)
        exportFieldToParaview(field_lidar,'lidar');
    end
    
    xmax = field_lidar.x0(1) + field_lidar.nx*field_lidar.hx;
    ymax = field_lidar.x0(2) + field_lidar.ny*field_lidar.hy;
    field_lidar.xf = [xmax,ymax];
end
    
if(~isempty(topoFieldName))
    disp('Read topo field')
    if(exportTopoParaview)
        exportFieldToParaview(field_topo,'topo');
    end
    
    xmax = field_topo.x0(1) + field_topo.nx*field_topo.hx;
    ymax = field_topo.x0(2) + field_topo.ny*field_topo.hy;
    field_topo.xf = [xmax,ymax];
else
    error('right now I need it, but I can copy the other')
end

tic; fprintf('Read height field...\n')
if(~isempty(lidarFieldName))
    %[field] = readHeighField(heighFieldName);
else
    [field_lidar] = readHeighField(mesh_facade);
end
elapsedTime = toc; fprintf('   ...%4.1f sec\n',elapsedTime);

field_lidar.x0 = field_lidar.x0 -translationVector;
field_lidar.xf = field_lidar.xf -translationVector;
field_topo.x0  = field_topo.x0  -translationVector;
field_topo.xf  = field_topo.xf  -translationVector;
translationVector_saved = translationVector; %%%Be careful: hack with transVect/transVect_saved
translationVector = [0.0 0.0];

fields.structured = true;
fields.ground = field_topo;
fields.roof = field_lidar;

