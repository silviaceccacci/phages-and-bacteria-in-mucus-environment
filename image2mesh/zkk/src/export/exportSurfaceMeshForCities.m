function []=exportSurfaceMeshForCities(surfaceMesh,fname)

    stream = fopen(['./Alya/' fname],'w');
    
%     fprintf(stream,'minX %f\n', min(mesh2D.X(:,1)));
%     fprintf(stream,'maxX %f\n', max(mesh2D.X(:,1)));
%     fprintf(stream,'minY %f\n', min(mesh2D.X(:,2)));
%     fprintf(stream,'maxY %f\n', max(mesh2D.X(:,2)));
%         
    fprintf(stream,'%d\n',size(surfaceMesh.X,1));
    fprintf(stream,'%d %f %f %f\n', [1:size(surfaceMesh.X,1); surfaceMesh.X']);
    
    fprintf(stream,'%d %d\n',size(surfaceMesh.T,1),size(surfaceMesh.T,2));
    fprintf(stream,'%d %d %d %d\n', [ 1:size(surfaceMesh.T,1); surfaceMesh.T']);
    
    fclose(stream);
    
end