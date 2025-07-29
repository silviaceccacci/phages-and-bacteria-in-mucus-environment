
function [] = printEdgesToPreserve(stream,mesh)

element.type = 'tri';
element.numVertices = 3; element.numEdges = 3; element.order = 1;
[boundNodes,bElems,bEdges] = giveBoundaryFromConnectivity(...
    mesh.T,size(mesh.X,2),element);


numEdges = length(bElems);
fprintf(stream,'\nEdges\n%d\n',numEdges);
for iedge=1:numEdges
    elemId      = bElems(iedge);
    localEdgeId = bEdges(iedge);
    nodes = getEdge(mesh.T(elemId,:),element,localEdgeId);
    fprintf(stream,'%d %d 0\n', nodes);
end
fprintf(stream,'\n');

fprintf(stream,'\nRequiredEdges %d\n',numEdges);
for iedge=1:numEdges
    fprintf(stream,'%d\n',iedge );
end
fprintf(stream,'\n');

end
