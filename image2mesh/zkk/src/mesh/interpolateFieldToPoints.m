function [values,elements,shapeF]=interpolateFieldToPoints(field,X)

    field.X = field.points;
    [elements,shapeF,typeCoord]=findPointsInMesh(field,X);
        
    values = shapeF(:,1).*field.z(field.T(elements,1)) +...
             shapeF(:,2).*field.z(field.T(elements,2)) +...
             shapeF(:,3).*field.z(field.T(elements,3)) ;

end