function []=exportFieldToParaview(field,fileName)

disp('BEGIN exporting lidar to paraview')

    minSampledTopo = field.x0;
    maxSampledTopo = field.x0+[field.hx*field.nx field.hy*field.ny];

    LX = maxSampledTopo(1)-minSampledTopo(1);
    LY = maxSampledTopo(2)-minSampledTopo(2);
    hhx = LX/(size(field.z,1)-1);
    hhy = LY/(size(field.z,2)-1);
    [XXX,YYY]=meshgrid(minSampledTopo(1):hhx:maxSampledTopo(1),minSampledTopo(2):hhy:maxSampledTopo(2));
    XXX = XXX';
    YYY = YYY';
   
    ZZZ = field.z;
    
    
    %tri = delaunay(XXX(:),YYY(:));
    %mesh.T = tri;
    
    npx = size(XXX,1);
    npy = size(XXX,2);
    numQuads = (npx-1)*(npy-1);
    T = zeros(numQuads,4);
    countElem = 0;
    for i=1:(npx-1)
        for j=1:(npy-1)
            countElem = countElem+1;
            n1 = i + (j-1)*npx;
            n2 = i+1 + (j-1)*npx;
            n3 = i+1 + j*npx;
            n4 = i + j*npx;
            T(countElem,:) = [ n1 n2 n3 n4];
        end
    end
    
    mesh.T = T;
    
    mesh.X = [XXX(:) YYY(:) ZZZ(:)];
    
    element.dim = 2;
    element.numCoord = 3;
    element.type = 'quad';
    element.order = 1;
    element.distribution = 'lineal';
    mesh.element = defineElement(element);
    
    mesh.name = fileName;
    %exportTriMeshToParaview(mesh);
    exportQuadMeshToParaview(mesh);
    
disp('END exporting lidar to paraview')
end