function [Z]=projectMeshToGround(mesh,field)

getElementsOverAverage = true;

dirichletBoundary = false;
minimumAtGround = false;

if(minimumAtGround && dirichletBoundary)
    error('one or the other')
end

%forceIntegration = 'gaussPoints';%'pixels','gaussPoints';
global fieldIntegrationType;
forceIntegration =fieldIntegrationType;

if(nargin<3)
    fakeBuildingHeight = 0.0;
end
%%
groundRegion = mesh.groundRegion;

if(isstruct(field))
    field_ground = field.ground;
    field_roof = field.roof;
else
    field_ground = field;
    field_roof = field;
end

%%
[shapeFunctions,gaussWeights,gaussPoints]=...
    computeShapeFunctionsReferenceElement(1,mesh.element.coord,mesh.element.orderQuadrature,mesh.element);

mesh.shapeFunctions = shapeFunctions(:,:,1);
mesh.gaussWeights = gaussWeights;
mesh.gaussPoints = gaussPoints; % we only need the functions, not the derivatives

numNodes = size(mesh.X,1);
Z = zeros(numNodes,1);

if(minimumAtGround)
    Z(:) = inf; %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Classify elements into roofs or ground
X2D = mesh.X(:,1:2);

groundElements = 1:size(mesh.T,1);
% groundElements = find(mesh.elementField==groundRegion);% from now on I'll be assuming that the interestElements are the first ones
groundMesh = mesh;
groundMesh.T = mesh.T(groundElements,:);

% roofElements = find(mesh.elementField>groundRegion);% from now on I'll be assuming that the interestElements are the first ones
% roofMesh = mesh;
% roofMesh.T = mesh.T(roofElements,:);


switch forceIntegration
    case 'pixels';
        disp('projecting...')
%         [projection] = project(field, interestMesh,'fieldToMesh'); % (ngauss,nelem)
%         [projection_roof] = project(field_roof, roofMesh,'fieldToMesh'); % (ngauss,nelem)
%         [projection_ground] = project(field_ground, groundMesh,'fieldToMesh'); % (ngauss,nelem)
%         disp('end projecting...')
% 
%         Zg = projection.z;
%         SFg = projection.shapeF;
%         elemsProj = projection.elements;
%         shapeFProj = projection.shapeF;
%               
%         shapeF_byElement = projection.shapeF_byElement;
%         numPointsInElem = shapeF_byElement.numPointsInElem;
%         shapeF_inElem = shapeF_byElement.shapeF;
%         points_inElem = shapeF_byElement.points;
        
    case 'gaussPoints'
        %disp('projectMeshToField_elemStruct:project only some elements to each field: roof to lidar, ground to topo')
        xg_g = mesh.X(groundMesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
               mesh.X(groundMesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
               mesh.X(groundMesh.T(:,3),1)*mesh.shapeFunctions(3,:);
        yg_g = mesh.X(groundMesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
               mesh.X(groundMesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
               mesh.X(groundMesh.T(:,3),2)*mesh.shapeFunctions(3,:);
        [Zg_ground] = project(field_ground, [xg_g(:) yg_g(:)],'pointsToField'); % (ngauss,nelem)
        Zg_ground = reshape(Zg_ground,size(xg_g));
        
%         xg_r = mesh.X(roofMesh.T(:,1),1)*mesh.shapeFunctions(1,:)+...
%                mesh.X(roofMesh.T(:,2),1)*mesh.shapeFunctions(2,:)+...
%                mesh.X(roofMesh.T(:,3),1)*mesh.shapeFunctions(3,:);
%         yg_r = mesh.X(roofMesh.T(:,1),2)*mesh.shapeFunctions(1,:)+...
%                mesh.X(roofMesh.T(:,2),2)*mesh.shapeFunctions(2,:)+...
%                mesh.X(roofMesh.T(:,3),2)*mesh.shapeFunctions(3,:);
%         [Zg_roof] = project(field_roof, [xg_r(:) yg_r(:)],'pointsToField'); % (ngauss,nelem)
%         Zg_roof = reshape(Zg_roof,size(xg_r));
         
        wg = gaussWeights/2.0;
  
    otherwise
        error('not implemented integration for the force term')
end

%% Least-squares minimization for the ground
projectionType_ground = 'minimization'; %minimization,interpolation

if(strcmp(projectionType_ground,'minimization'))
    %groundElements = interestElements(find(mesh.elementField(interestElements)==groundRegion));
    im = zeros(length(groundElements)*size(mesh.T,2)^2,1);
    jm = zeros(length(groundElements)*size(mesh.T,2)^2,1);
    km = zeros(length(groundElements)*size(mesh.T,2)^2,1);
    mloc1 = [ 1 2 3 1 2 3 1 2 3];
    mloc2 = [ 1 1 1 2 2 2 3 3 3];
    count = 0 ;
    F = zeros(numNodes,1);
    groundNodes = zeros(numNodes,1);
    for locelem = 1:length(groundElements)
        ielem = groundElements(locelem);
        
        elemNodes = mesh.T(ielem,:);
        groundNodes(elemNodes) = 1;
        elemArea   = computeArea(mesh.X(elemNodes,1:2));

        Mdiag = (mesh.shapeFunctions.*mesh.shapeFunctions)*mesh.gaussWeights;
        Melem = diag(Mdiag);
        Melem(1,2) = (mesh.shapeFunctions(1,:).*mesh.shapeFunctions(2,:))*mesh.gaussWeights;
        Melem(1,3) = (mesh.shapeFunctions(1,:).*mesh.shapeFunctions(3,:))*mesh.gaussWeights;
        Melem(2,3) = (mesh.shapeFunctions(2,:).*mesh.shapeFunctions(3,:))*mesh.gaussWeights;
        Melem(3,2) = Melem(2,3);
        Melem(3,1) = Melem(1,3);
        Melem(2,1) = Melem(1,2);
        Melem = Melem * elemArea;

        if(strcmp(forceIntegration,'pixels') && numPointsInElem(ielem)==0)
            mesh.X(mesh.T(ielem,:),:)
            elemArea
            figure; grid on; hold on;
            plot(mesh.X(mesh.T(ielem,[1 2 3 1]),1),mesh.X(mesh.T(ielem,[1 2 3 1]),2))
            plot(field.points(:,1),field.points(:,2),'*')
            error('No pixels in this element????')
        end

        
        switch forceIntegration
            case 'pixels'        
                sfaux = shapeF_inElem{ielem};
                zaux = Zg( points_inElem{ielem} );

                Felem = bsxfun(@times,sfaux,zaux);
                Felem = sum(Felem,1)'/numPointsInElem(ielem)/2.0;
                Felem = Felem * elemArea ;       
            case 'gaussPoints'
                Felem = bsxfun(@times,mesh.shapeFunctions,Zg_ground(locelem,:))*mesh.gaussWeights/2.0;
                Felem = Felem * elemArea;
                F(elemNodes) = F(elemNodes) + Felem;
                
                if(minimumAtGround)
                    minZ_elem = min(Zg_ground(locelem,:));
                    for iinode = 1:length(elemNodes)
                       if( Z(elemNodes(iinode)) >minZ_elem ) 
                           Z(elemNodes(iinode)) = minZ_elem;
                       end
                    end
                end
        end
                 
        F(elemNodes) = F(elemNodes) + Felem;

        im((count+1):(count+length(elemNodes)^2)) = elemNodes(mloc1);
        jm((count+1):(count+length(elemNodes)^2)) = elemNodes(mloc2);
        km((count+1):(count+length(elemNodes)^2)) = Melem(:);
        count = count + length(elemNodes)^2;

        % EI!!!! PROVA!!!
        if(dirichletBoundary)
            Z(elemNodes) = min(Zg(ielem,:)); %%%%%%%%%%%%%%% EIEIEIEIEIEIEIE
        end
    end
    groundNodes = find(groundNodes);
    
    if(dirichletBoundary)
        disp('Fixing boundary (dirichlet)')
        boundaryNodes = [mesh.boundaryNodes.exterior ; mesh.boundaryNodes.ground];%%%%%%%%%%%%%%% EIEIEIEIEIEIEIE
        groundNodes = setdiff(groundNodes,boundaryNodes);%%%%%%%%%%%%%%% EIEIEIEIEIEIEIE
    end
    
    M = sparse(im,jm,km);
    
    M = M(groundNodes,:);
    
    if(dirichletBoundary)
        Fsubstract = M(:,boundaryNodes)*Z(boundaryNodes);%%%%%%%%%%%%%%% EIEIEIEIEIEIEIE
    end
    
    M = M(:,groundNodes);
    
    F = F(groundNodes);
    
    if(dirichletBoundary)
        F = F-Fsubstract;%%%%%%%%%%%%%%% EIEIEIEIEIEIEIE
    end
    
    if(~minimumAtGround)
        fprintf('      Condition number of the matrix: %1.2e\n',condest(M))
        Z(groundNodes) = M\F;
    end
    
    %%Z(groundNodes) = F;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
elseif(strcmp(projectionType_ground,'interpolation'))
    error('dont use this')
    groundNodes = zeros(numNodes,1);
    nodeNeighElem = zeros(size(mesh.X,1),1);
    %for ielem = 1:size(mesh.T,1)
    for locelem = 1:length(interestElements)
        ielem = interestElements(locelem);
        iregion = mesh.elementField(ielem);
        if(iregion==groundRegion)
            elemNodes = mesh.T(ielem,:);
            elemHeight = Zg(ielem,:)*wg;
            Z(elemNodes) = Z(elemNodes) + elemHeight;
            nodeNeighElem(elemNodes) = nodeNeighElem(elemNodes) + 1;

            groundNodes(elemNodes) = 1;
        end
    end
    groundNodes = find(groundNodes);
    Z(groundNodes) = Z(groundNodes)./nodeNeighElem(groundNodes);
         
end
%[min(Z(groundNodes)) max(Z(groundNodes))]


%% output
% X3D = [X2D,  Z];
% 
% mesh_facade3D = mesh;
% mesh_facade3D.X = X3D;
% mesh_facade3D.name = [mesh.name '3D'];
% 
% mesh_facade3D.holeBuildings =[];
% %mesh_facade3D.holeBuildings = holeBuildings;
% 
% % figure(234)
% % listP = 1:500:size(field.points,1);
% % plot3(field.points(listP,1),field.points(listP,2),field.z(listP),'*')
% % hold on
% % plot3(X3D(:,1),X3D(:,2),X3D(:,3),'ro')
end

function [area] = computeArea(X)

if(size(X,1)==3)
    v1 = [ X(2,:)-X(1,:) , 0];
    v2 = [ X(3,:)-X(1,:) , 0];
    area = norm(cross(v1,v2))/2.0;
else
    error('not implemented')   
end

end