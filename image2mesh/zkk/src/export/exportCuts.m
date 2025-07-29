function [meshTopo]=exportCuts(fileName,field,mesh2D,mesh_facade3D,writeFile,translationVector_saved)

%     if(field.nx<=2) % if there is no topography
%         
%         disp('No DEM to export, h=0')
%         
%         %meshTopo = mesh2D;
%         
%         myEps = 1;
%         x0 = min(mesh2D.X(:,1)) - myEps;
%         y0 = min(mesh2D.X(:,2)) - myEps;
%         x1 = max(mesh2D.X(:,1)) + myEps;
%         y1 = max(mesh2D.X(:,2)) + myEps;
%         
%         X = [x0 y0; x1 y0 ; x1 y1; x0 y1];
%         T = [ 1 2 4; 2 3 4];
%         meshTopo.X = [X zeros(size(X,1),1)];
%         meshTopo.T = T;
%         
%     else
        z = zeros(size(mesh2D.X,1),1);
        groundElements = find(mesh_facade3D.elementField==mesh_facade3D.groundRegion);
        groundNodes = unique(mesh_facade3D.T(groundElements,:));
        z(groundNodes)=mesh_facade3D.X(groundNodes,3);
        notAssignedNodes = setdiff(1:size(mesh2D.X,1),groundNodes);
        % HE CANVIAT AIXO XQ NO TRASLLADAVA EL PUNT: 25/03/19
        %z(notAssignedNodes) = project(field, mesh2D.X(notAssignedNodes,:),'pointsToField');
        z(notAssignedNodes) = project(field, mesh2D.X(notAssignedNodes,:)+translationVector_saved,'pointsToField');

        meshTopo.X = [mesh2D.X z];
        meshTopo.T = mesh2D.T;
%     end
    
    if(writeFile)
        exportSurfaceMeshForCities(meshTopo,[fileName '.cutSurface']);
    end
    
    
%     figure(12121)
%     plot3(mesh2D.X(:,1),mesh2D.X(:,2),z,'.')
%     axis equal

% global cmglob;
% mesh2D.X = [mesh2D.X z];
% mesh2D.X = bsxfun(@minus,mesh2D.X,cmglob);
% mesh2D.name = 'cutSrf';
% mesh2D.element = mesh_facade3D.element;
% exportTriMeshToParaview(mesh2D);


end

%     %mesh2D.elementField = mesh2D.groundRegion;
% 
%     %z = project(field, mesh2D.X,'pointsToField');
%     z = zeros(size(mesh2D.X,1),1);
%     groundElements = find(mesh_facade3D.elementField==mesh_facade3D.groundRegion);
%     groundNodes = unique(mesh_facade3D.T(groundElements,:));
% %     groundNodes = (mesh_facade3D.T(groundElements,:));
% %     groundNodes = groundNodes(:);
%     z(groundNodes)=mesh_facade3D.X(groundNodes,3);
%     notAssignedNodes = setdiff(1:size(mesh2D.X,1),groundNodes);
%     z(notAssignedNodes) = project(field, mesh2D.X(notAssignedNodes,:),'pointsToField');
%     
%     meshTopo.X = [mesh2D.X z];
%     meshTopo.T = mesh2D.T;
%     
%     if(writeFile)
%         exportSurfaceMeshForCities(meshTopo,[fileName '.cutSurface']);
%     end
%     
%     
% %     figure(12121)
% %     plot3(mesh2D.X(:,1),mesh2D.X(:,2),z,'.')
% %     axis equal
% 
% % global cmglob;
% % mesh2D.X = [mesh2D.X z];
% % mesh2D.X = bsxfun(@minus,mesh2D.X,cmglob);
% % mesh2D.name = 'cutSrf';
% % mesh2D.element = mesh_facade3D.element;
% % exportTriMeshToParaview(mesh2D);
