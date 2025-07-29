function [polylines,translationVector]=read1DGeometry(fileName,geomType,fileName_point)

if(strcmp(geomType,'gid'))
    [polylines]=read1DGeometry_gid(fileName);
else
    error('Not know geometry type')
end

% if(~isempty(fileName_point))
%     stream = fopen(fileName_point,'r');     
%     if (stream ~= -1)  
%         
%     end
% end
  
    fid = fopen([fileName_point '.txt'], 'r');    
    translationVector = fscanf(fid, '%f',2)';
    fclose(fid);

%    fixTrans = 1.0e+03 *[1.198559451415902  -1.284981071637943];
%    translationVector=translationVector + fixTrans;
end

function [polylines]=read1DGeometry_gid(fileName)

    fid = fopen(fileName, 'r');
    
    fgetl(fid);
    fgetl(fid);
    
    X = zeros(3,0);
    newNode = 0;
    word = fscanf(fid,'%s',1);
    while(strcmp(word,'End')==false)
        newNode = newNode+1;
        X(1:3,newNode) = fscanf(fid, '%f %f %f\n',3);
        word = fscanf(fid,'%s',1);
%         if(mod(newNode,10000)==0)
%             fprintf('%d ',newNode)
%         end
    end
%     fprintf('\n ')
    X = X(1:2,:);
    fgetl(fid);
    fgetl(fid);
    fgetl(fid);

    T = zeros(2,0);
    newElem = 0;
    word = fscanf(fid,'%s',1);
    while(strcmp(word,'End')==false)
        newElem = newElem+1;
        T(1:2,newElem) = fscanf(fid, '%d %d\n',2);
        word = fscanf(fid,'%s',1);
%         if(mod(newElem,10000)==0)
%             fprintf('%d ',newElem)
%         end
    end
    fclose(fid);

    polylines.T = T';
    polylines.X = X';
    
end