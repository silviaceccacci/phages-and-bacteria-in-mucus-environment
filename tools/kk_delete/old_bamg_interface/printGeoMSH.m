function printGeoMSH(fullName,mesh,h)

    stream = fopen(fullName,'w');
    
    fprintf(stream,'MeshVersionFormatted 0\n\n');
    fprintf(stream,'Dimension 2\n');
    fprintf(stream,'\nVertices\n%d\n',size(mesh.X,2));
    fprintf(stream,'%f %f %d\n', [mesh.X; zeros(1,size(mesh.X,2))]);
    fprintf(stream,'\nEdges\n%d\n',size(mesh.T,1));
    fprintf(stream,'%d %d %d\n', [ mesh.T';  zeros(1,size(mesh.T,1))]);
    
    if(nargin>2)
        fprintf(stream,'\nhVertices\n');
        fprintf(stream,'%f ', h);
        fprintf(stream,'\n');
    end
    
    fprintf(stream,'\nEnd\n');
    
    fclose(stream);

end