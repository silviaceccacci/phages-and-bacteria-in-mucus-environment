function [z_on_x] = project_toMesh(field, mesh, options)  
        
    switch options
        case 'nodes'
            [z_on_x] = project_toMeshNodes(field, mesh);
        case 'gauss'
            [z_on_x] = project_toMeshGauss(field, mesh);
        otherwise
            error('not implemented')
            
    end

    % list the gauss points xg on all the physical elements in the ground or roofs
    % [z_on_xg] = interpolateFromMesh(field, xg)
    % compute: 
    %   - case roof: integral on all the region if the height to get const
    %   height
    %   - case floor: use to get least squares minimizing surface

    %error('Hey I am here and I will kill this in five minuts! ;-p!')
    
end

function [z_on_x] = project_toMeshNodes(field, mesh)
    z_on_x = griddata(field.points(:,1),field.points(:,2),field.z(:),mesh.X(:,1),mesh.X(:,2));
end

function [z] = project_toMeshGauss(field, mesh)   

    %x = (mesh.X(mesh.T(:,1),1)+mesh.X(mesh.T(:,2),1)+mesh.X(mesh.T(:,3),1))/3.0;
    %y = (mesh.X(mesh.T(:,1),2)+mesh.X(mesh.T(:,2),2)+mesh.X(mesh.T(:,3),2))/3.0;

    x = mesh.X(mesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
        mesh.X(mesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
        mesh.X(mesh.T(:,3),1)*mesh.shapeFunctions(3,:);
    y = mesh.X(mesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
        mesh.X(mesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
        mesh.X(mesh.T(:,3),2)*mesh.shapeFunctions(3,:);
    
    x = x(:);
    y = y(:);
    
%     figure(100)
%     hold on
%     plot(mesh.X(:,1),mesh.X(:,2),'o')
%     plot(x,y,'*')
%     pause()
    
    %wg = 1;%ones(size(mesh.T,1),1);
    wg = mesh.gaussWeights;

    %z.z = griddata(field.points(:,1),field.points(:,2),field.z(:),x,y);
    [values,z.elements,z.shapeF] =interpolateFieldToPoints(field,[x y]);       
    
    z.z = reshape(values,size(mesh.T,1),size(mesh.shapeFunctions,2));
    
    z.wg = wg; % because these shapeF are in [0,1][0,1]
end