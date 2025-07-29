
function []=exportAlya(mesh,fileName,options)

    write_dom(mesh,fileName);
    write_geo(mesh,fileName);
    write_fix(mesh,fileName);
    write_BCS(mesh,fileName,options);

end

function [] = write_dom(mesh,fileName)             
fname = [fileName '.dom.dat'];
stream = fopen(['./Alya/' fname],'w');
  %
out_string = [...
   '$------------------------------------------------------------\n' ...
   'WIND_MESH\n' ...
   '  SURFACE_POINTS %d\n' ...
   '  NUMBER_OF_CELLS %d\n' ...
   '  FIRST_CELL_SIZE %f\n' ...
   'END_WIND_MESH\n' ...
   '$------------------------------------------------------------\n' ...
   'DIMENSIONS\n' ...
   '  NODAL_POINTS  %d\n' ...
   '  ELEMENTS      %d\n' ...
   '  SPACE_DIMENSIONS   3\n' ...
   '  TYPES_OF_ELEMENTS  TET04\n' ...
   '  BOUNDARIES    %d\n' ...
   'END_DIMENSIONS\n' ...
   '$------------------------------------------------------------\n' ...
   'STRATEGY\n' ...
   '  INTEGRATION_RULE           Closed\n' ...  
   '  DOMAIN_INTEGRATION_POINTS  8\n' ...
   'END_STRATEGY\n' ...
   '$-------------------------------------------------------------\n' ...
   'GEOMETRY\n' ...
   '  GROUPS = 500\n' ...
   '  TYPES, ALL=TET04\n' ...
   '  END_TYPES\n' ...
   '  INCLUDE %s\n' ...
   '  INCLUDE %s\n' ...
   'END_GEOMETRY\n' ...  
   '$-------------------------------------------------------------\n' ...
   'SETS\n' ...
   'END_SETS\n' ...
   '$-------------------------------------------------------------\n' ...
   'BOUNDARY_CONDITIONS\n' ...
   '  INCLUDE %s\n' ...  
   '  GEOMETRICAL_CONDITIONS\n' ...
   '     FREESTREAM   1,2,3,4,6\n' ...
   '$     SYMMETRY     6\n' ...
   '     WALL_LAW     5\n' ...
   '$     ANGLE, GEO_ANGLE = 150.0\n' ...
   '     CRITERION_FREESTREAM : VALUE_FUNCTION =  1 toler 90.0\n' ...
   '  END_GEOMETRICAL_CONDITIONS\n' ...
   '  VALUES, FUNCTION=1, DIMENSION=3\n' ...
   '     INCLUDE %s\n' ...
   '  END_VALUES\n' ...
   '  VALUES, FUNCTION=2, DIMENSION=1\n' ...
   '     INCLUDE %s\n' ...
   '  END_VALUES\n' ...
   '  VALUES, FUNCTION=3, DIMENSION=1\n' ...
   '     INCLUDE %s\n' ...
   '  END_VALUES\n' ...
   'END_BOUNDARY_CONDITIONS\n' ...
   '$-------------------------------------------------------------'];
fprintf(stream,out_string,...
    0,0,0,... %npoin2d,ncellz,dz0, ...
    size(mesh.X,1),size(mesh.T,1),size(mesh.boundaries.faces,1), ... 
    [fileName '.geo'], ...
    [fileName '.zo.dat'], ...
    [fileName '.fix'], ...
    ['./AlyaFix/' fileName '.nsi.ini'], ...   
    ['./AlyaFix/' fileName '.tur1.ini'], ...   
    ['./AlyaFix/' fileName '.tur2.ini']);
fclose(stream);
end

function [] = write_geo(mesh,fileName)  
    %
    %*** 1. Geometry (mesh) file 
    %
    fname = [fileName '.geo'];
    stream = fopen(['./Alya/' fname],'w');
    %
    %*** Writes nodal connectivities
    %
    fprintf(stream,'ELEMENTS\n');
    fprintf(stream,'%d %d %d %d %d\n', [ 1:size(mesh.T,1); mesh.T']); %(:,[1 3 2 4])
    fprintf(stream,'END_ELEMENTS\n');
    %
    %*** Writes coordinates
    %
    fprintf(stream,'COORDINATES\n');
    fprintf(stream,'%d %f %f %f\n', [1:size(mesh.X,1); mesh.X']);
    fprintf(stream,'END_COORDINATES\n');
    %
    %*** Writes boundaries nodal connectivities
    %
    fprintf(stream,'BOUNDARIES, ELEMENTS\n');%---> Note the order, opposite to OF. Rotation inwards
    fprintf(stream,'%d %d %d %d %d\n', ...
        [ 1:length(mesh.boundaries.marks); mesh.boundaries.faces(:,[1 3 2])' ; mesh.boundaries.elems']);
    fprintf(stream,'END_BOUNDARIES\n');
    %
    %*** Writes skew systems
    %
    fprintf(stream,'SKEW_SYSTEMS\n');
    fprintf(stream,'END_SKEW_SYSTEMS\n');
    %
    fclose(stream);
    %
  
end

function [] = write_fix(mesh,fileName)

    fname = [fileName '.fix'];
    stream = fopen(['./Alya/' fname],'w');
    
    fprintf(stream,'ON_BOUNDARIES\n');
    fprintf(stream,'%d %d\n', ...
        [ 1:length(mesh.boundaries.marks); mesh.boundaries.marks']);
    fprintf(stream,'END_ON_BOUNDARIES\n');
    
    fclose(stream);
    
end

function [] = write_BCS(mesh,fileName,options)
    %% Compute nodes height for initial condition
    i = 1;
    projection = [];
    %projection{i} = 'dem'; i = i+1;
    %projection{i} = 'lidar'; i = i+1;
    %projection{i} = 'cityMesh'; i = i+1;
    %projection{i} = 'topoMesh';
    
    floorNodes = find(mesh.boundaries.marks==5);
    
    for i=1:length(projection)
        switch projection{i}
            case 'dem'
                
                field_topo = mesh.fields.ground;
                zProj = project(field_topo,mesh.X(:,1:2),'pointsToField');

                mesh.heightTopo = mesh.X(:,3)-zProj;
                
                negativeNodes = find(mesh.heightTopo<=0.0);
                mesh.heightTopo(negativeNodes) = 0.0;
                
            case 'lidar'
                field_lidar = mesh.fields.roof;
                field_topo = mesh.fields.ground;
                xnodes = mesh.X(:,1);
                ynodes = mesh.X(:,2);
                switch mesh.domainLimits.type
                    case {'intersection','cadastre','lidar'}
                        zProj = project(field_lidar,mesh.X(:,1:2),'pointsToField');
                    otherwise % we need to find which nodes to project to lidar and which to dem
                        pointsToLidar = 1:size(mesh.X,1);
                        pointsToLidar = pointsToLidar(find(xnodes(pointsToLidar)>field_lidar.x0(1)));
                        pointsToLidar = pointsToLidar(find(xnodes(pointsToLidar)<field_lidar.xf(1)));
                        pointsToLidar = pointsToLidar(find(ynodes(pointsToLidar)>field_lidar.x0(2)));
                        pointsToLidar = pointsToLidar(find(ynodes(pointsToLidar)<field_lidar.xf(2)));
                        pointsToDem = setdiff(1:size(mesh.X,1),pointsToLidar);
                        zProj = zeros(size(mesh.X,1),1);
                        zProj(pointsToLidar) = project(field_lidar,mesh.X(pointsToLidar,1:2),'pointsToField');
                        if(~isempty(pointsToDem))
                            zProj(pointsToDem)   = project(field_topo,mesh.X(pointsToDem,1:2),'pointsToField');
                        end
                end

                mesh.heightCity = mesh.X(:,3)-zProj;

                mesh.heightCity(floorNodes) = 0.0;
            
                negativeNodes = find(mesh.heightCity<=0.0);
                mesh.heightCity(negativeNodes) = 0.0;
                
            case 'topoMesh'
                
                fprintf('Computing heights with respect to topo mesh\n')
                
                optionsProj = 'nodes';%'gauss'
                
                clear field;
                field.structured = false;
                field.points = options.meshTopo.X(:,1:2);
                field.z = options.meshTopo.X(:,3);
                heightGround = project(field, mesh, optionsProj);
                
                mesh.heightTopo = mesh.X(:,3)-heightGround;
                
                mesh.heightTopo(floorNodes) = 0.0;
                
                negativeNodes = find(mesh.heightTopo<0.0);
                mesh.heightTopo(negativeNodes) = 0.0;
                
                %error('not programmed')
                
            case 'cityMesh'
                
                fprintf('Computing heights with respect to city mesh\n')
                
                optionsProj = 'nodes';%'gauss'
                
                clear field;
                field.structured = false;
                x = options.meshFacades.X(:,1:2);
                z = options.meshFacades.X(:,3);
%                 heightGround = project(field, mesh, optionsProj);
                heightGround = griddata(x(:,1),x(:,2),z(:),mesh.X(:,1),mesh.X(:,2));%,'nearest');


%                 field.X = options.meshFacades.X;
%                 field.T = options.meshFacades.T;
%                 optionsFindPoints.bin = true;
%                 optionsFindPoints.numBins = max(1,round(sqrt(size(field.T,1))/100));
%                 [elements,shapeF,typeCoord]=findPointsInMesh(field,mesh.X(:,1:2),optionsFindPoints);
% 
%                 heightGround= shapeF(:,1).*field.z(field.T(elements,1)) +...
%                               shapeF(:,2).*field.z(field.T(elements,2)) +...
%                               shapeF(:,3).*field.z(field.T(elements,3)) ;
%                 
                mesh.heightCity = mesh.X(:,3)-heightGround;
                
                mesh.heightCity(floorNodes) = 0.0;
                
                negativeNodes = find(mesh.heightCity<0.0);
                mesh.heightCity(negativeNodes) = 0.0;
                %error('not programmed')
                
        end
    end


    %% Export
    fname = [fileName '.bcs'];
    stream = fopen(['./Alya/' fname],'w');
    
    numNodes = size(mesh.X,1);
    if(isfield(mesh,'rough'))
        rough = mesh.rough;
    else
        rough = ones(numNodes,1);
    end
%     if(isfield(mesh,'height'))
%         height = mesh.height;
%     else
%         height = zeros(numNodes,1);
%     end
    
    fprintf(stream,'%d\n',numNodes);    
    if(length(projection)==0)
        fprintf(stream,'%f\n', ...
            rough');
        
        exportSurfaceMeshForCities(options.meshTopo,[fileName '.surfaceDEM']);
                                
%         meshTopoForUs = options.meshTopo;
%         meshTopoForUs.name = 'mallaDEM';
%         exportTriMeshToParaview(meshTopoForUs)
                
        roofElements = (options.meshFacades.elementField>options.meshFacades.groundRegion);
        options.meshFacades.T = options.meshFacades.T(roofElements,:);
        exportSurfaceMeshForCities(options.meshFacades,[fileName '.surfaceRoof']);
        
    elseif(length(projection)==1)
        fprintf(stream,'%f %f %f\n', ...
            [ rough'; mesh.height']);
    elseif(length(projection)==2)
        fprintf(stream,'%f %f %f\n', ...
            [ rough'; mesh.heightCity'; mesh.heightTopo']);        
    else
        error('Something wrong')
    end
    fclose(stream);
    
end




