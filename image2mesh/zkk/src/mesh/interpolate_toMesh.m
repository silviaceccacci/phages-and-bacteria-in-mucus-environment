function [z_on_x] = interpolate_toMesh(field, x)  
    error('should not be used any more')
    %% Node interpolation (fails on discontinuities)
    %z_on_x = griddata(field.pixelPoints(:,1),field.pixelPoints(:,2),field.z(:),x(:,1),x(:,2));
    
    %% Interpolation on the center of the triangle
    z_on_xg = pdeintrp(field.inventedField.pixelPoints,field.T,field.z);
    
    %% Interpolate on mesh: does not work
%     p = [field.pixelPoints(:,1),field.pixelPoints(:,2)]';
%     t = field.T';%[ field.T ones(size(field.T,1),1)]';
%     z_on_x = tri2grid(p,t,field.z(:),x(:,1),x(:,2));
    
    %% hand made


end