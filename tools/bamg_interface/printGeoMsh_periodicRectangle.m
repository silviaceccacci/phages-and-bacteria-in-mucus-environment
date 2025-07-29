function printGeoMsh_periodicRectangle(fullName,Xcorners)

%        3
%   n4 ----- n3
%   |        |
% 4 |        | 2
%   |        |
%   n1 ----- n2
%        1


    stream = fopen(fullName,'w');
    
    fprintf(stream,'MeshVersionFormatted 0\n\n');
    fprintf(stream,'Dimension 2\n');
    
    fprintf(stream,'\nVertices\n%d\n',size(Xcorners,2));
    fprintf(stream,'%f %f %d\n', [Xcorners; zeros(1,size(Xcorners,2))]);

    fprintf(stream,'\nEdges\n4\n');
    fprintf(stream,'1 2 0\n');
    fprintf(stream,'2 3 0\n');
    fprintf(stream,'3 4 0\n');
    fprintf(stream,'4 1 0\n');
    
% • EquivalencedEdges (I) NbOfEquivalencedEdges
% ( @@Edge1
% i
% , @@Edge2
% i
% , i=1 , NbOfEquivalencedEdges)

    fprintf(stream,'\nEquivalencedEdges\n2\n');
    fprintf(stream,'1 3\n');
    fprintf(stream,'2 4\n');

    fprintf(stream,'\nEnd\n');
    
    fclose(stream);

end