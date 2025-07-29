function [mesh]=classifyElementRegionFromSTL(mesh,tolHeight)

tolHeightSTL = min(mesh.X(:,3))+tolHeight;%max(8,tolHeight);

numElements = size(mesh.T,1);

groundRegion = 1;
facadeRegion = 0;
roofRegion   = 2; % If cadastre, one roof region per roof. With STL, just one region for all.

tolZeroNormals = 5e-1;% 1e-1;%1e-2;
[normalElements]=computeNormalsMeshElements(mesh);
facadeElements = find( abs(normalElements(:,3))<tolZeroNormals );
nonFacadeElements = setdiff( 1:numElements,facadeElements );
roofElements = nonFacadeElements(find( mesh.X(mesh.T(nonFacadeElements,1),3)>tolHeightSTL ));

elementField = groundRegion*ones(numElements,1);
elementField(facadeElements) = facadeRegion;
elementField(roofElements) = roofRegion;

mesh.elementField = elementField;

mesh.groundRegion = groundRegion;
mesh.facadeRegion = facadeRegion;
mesh.roofRegion = roofRegion;

end

function [normalElements]=computeNormalsMeshElements(mesh)
    
    v1 = mesh.X(mesh.T(:,2),:) - mesh.X(mesh.T(:,1),:);
    v2 = mesh.X(mesh.T(:,3),:) - mesh.X(mesh.T(:,1),:);

    normalElements = [...
        v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2) , ...
        v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3) , ...
        v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1) ];
    
    normNormals = sqrt(sum( normalElements.^2 ,2));
    normalElements = bsxfun(@rdivide,normalElements,normNormals);

end

% facadeRegion = mesh.groundRegion - 1;
% newElementRegions = facadeRegion*ones(size(Tfacades,1),1);
% mesh.facadeRegion = facadeRegion;
% 
% mesh.facadeElementsRegion = facadeElementsRegion;
% 
% mesh.T = newT;
% mesh.X = newX;
% mesh.elementField = [mesh.elementField ; newElementRegions];
