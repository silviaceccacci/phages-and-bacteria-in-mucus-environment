function []=exportMetaData(fileName,mesh2D,cm)


    fname = [fileName '.metaData'];
    stream = fopen(['./Alya/' fname],'w');
    
    fprintf(stream,'x0 %f\n',   cm(1));
    fprintf(stream,'y0 %f\n',   cm(2));
    minX = min(mesh2D.X(:,1));
    fprintf(stream,'minX %f %f\n', minX, minX+cm(1));
    maxX = max(mesh2D.X(:,1));
    fprintf(stream,'maxX %f %f\n', maxX, maxX+cm(1));
    minY = min(mesh2D.X(:,2));
    fprintf(stream,'minY %f %f\n', minY, minY+cm(2));
    maxY = max(mesh2D.X(:,2));
    fprintf(stream,'maxY %f %f\n', maxY, maxY+cm(2));
        
    fclose(stream);
    

end