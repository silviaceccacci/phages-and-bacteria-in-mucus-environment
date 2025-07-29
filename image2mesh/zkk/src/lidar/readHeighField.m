function [field] = readHeighField(heighFieldName)
    if(isstruct(heighFieldName))
        [field] = generateMockField(heighFieldName);
    else
        [field] = readLidarToField(heighFieldName);
    end
end


function [field] = generateMockField(mesh_facade)

    % GENERATE MOCK IMAGE
    %lx = max(mesh_facade.X(:,1))-min(mesh_facade.X(:,1));
    %hx = lx/500;
    hx = 1;%2.5;%2.5;
    hy = hx;
    margin = 1*max([hx hy]);
    minX = min(mesh_facade.X(:,1))-margin;
    maxX = max(mesh_facade.X(:,1))+margin;
    minY = min(mesh_facade.X(:,2))-margin;
    maxY = max(mesh_facade.X(:,2))+margin;
    px = minX:hx:maxX;
    py = minY:hy:maxY;
    [Xg,Yg]=meshgrid(px,py);
    Xg=Xg';% per compensar el project al reves
    Yg=Yg';% per compensar el project al reves
    %Zg = (Xg-minX).*(Yg-minY)/((maxX-minX)*(maxY-minY));
    Zg = 400*((Xg-minX).^2 + (Yg-minY).^2 )/((maxX-minX).^2 + (maxY-minY).^2); %Xg/(maxX-minX) + Yg/(maxY-minY); %
    %Zg = sqrt(Xg.*Yg) ./sqrt(Xg.*Xg + Yg.*Yg);
%     Zg = (20.*Zg).^2;
%     %Zg = 20*sin(100*Zg/(maxX-minX)) + 20*Zg.^2;

    x0 = [minX minY];
    width = maxX-minX;
    height = maxY-minY;

    field.z = Zg;
    field.hx = hx;
    field.hy = hy;
    field.x0 = x0;
    %field.width = width;
    %field.height = height;
    %field.x = px;
    %field.y = py;
    field.points = [Xg(:) Yg(:)];

    %figure(10)
    %surf(Xg,Yg,Zg)
    
    
    field.structured = true;

end