function [projection] = project_structured(field, target, options)  

    switch options
        case 'pointsToField'
            [projection] = project_structured_pointsToField(field, target);
        case 'fieldToMesh'
            [projection] = project_structured_fieldToMesh(field, target);
        otherwise
            error('not implemented')
            
    end

%     [location] = locate(x, field);
%     
%     z = field.z;
%     [m, n, numComponents] = size(z);
%     
%     numPoints = size(x, 1);
%     
%     z_on_x = zeros(numPoints, numComponents);
%     
%     for c = 1:numComponents
%         zc = z(:,:,c);     
%         z_on_x(:, c) = zc(location.l);
%     end
    
end

function [z_on_x] = project_structured_pointsToField(field, x)  

    %disp('locate changed: return to original');
    [location] = locate(x, field);
    
    z = field.z;
    [m, n, numComponents] = size(z);
    
    numPoints = size(x, 1);
    
    z_on_x = zeros(numPoints, numComponents);
    
    for c = 1:numComponents
        zc = z(:,:,c);    
        z_on_x(:, c) = zc(location.l); % THIS IS THE CORRECT LINE!!!!!!!!!!!!!!!!!
        %z_on_x(:, c)  = 0.0;
        %myIndex = find(location.l>0);
        %myIndex = myIndex(find(location.l<=numel(zc)));
        %z_on_x(myIndex, c) = zc(location.l(myIndex));
    end
    
end

function [projection] = project_structured_fieldToMesh(field, mesh)  

    options.classifyByElement = true;
    options.stopWithNotFound = false;
    options.bin = true;
    options.numBins = 10;

    [elements,shapeF,typeCoord,shapeF_byElement]=findPointsInMesh(mesh,field.points,options);
    
    projection.elements = elements;
    projection.shapeF = shapeF;%(:,:,1);
    projection.z = field.z(:);
    projection.shapeF_byElement = shapeF_byElement;
    
end



