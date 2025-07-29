
function [] = printNodeAttributes_paraview(stream,meshLin,nodalFields)

    nOfFields = size(nodalFields,2);
    
    intVect = [nOfFields  ones(1,nOfFields)];
    
    fprintf(stream,'%d ',intVect(1:(length(intVect)-1)));
    fprintf(stream,'%d',intVect(length(intVect)));
    fprintf(stream,'\n');
    fprintf(stream,'qIni, None\n');
    fprintf(stream,'qGeo, None\n');
    fprintf(stream,'%d %f %f \n',[1:size(meshLin.X,2); nodalFields']);
    
end