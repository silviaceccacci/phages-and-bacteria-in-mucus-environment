
function [mesh2D] = getMesh2DFrom3DMesh(mesh)
    nelems = size(mesh.T,1);
    element = mesh.element;
    numNodesFace = giveNumNodesFaceFromOrder(element.order, element.type);
    T2D = zeros(nelems*element.numFaces, numNodesFace);
    qualityFaces = zeros(size(T2D,1), 1);
    if(isfield(mesh,'quality')==false)
        mesh.quality = 0.5*ones(size(mesh.T,1),1);
    end
    if(isfield(mesh,'qualities'))
       namesQ = fieldnames(mesh.qualities);
       for i=1:length(namesQ)
           qualitiesSubelements.(namesQ{i}) = zeros(size(qualityFaces));
       end
    end
    if(isfield(mesh,'distoritons'))
       namesQ = fieldnames(mesh.distoritons);
       for i=1:length(namesQ)
           distoritonsSubelements.(namesQ{i}) = zeros(size(qualityFaces));
       end
    end
    if(isfield(mesh,'errorDisc'))
       namesE = fieldnames(mesh.errorDisc);
       for i=1:length(namesE)
           errorDisc.(namesE{i}) = zeros(size(qualityFaces));
       end
    end

    CM = zeros(3,size(T2D,1));
    count = 0;
    for iElem = 1:nelems
        xyz = sum(mesh.X(:,mesh.T(iElem,:)),2)/element.numNod;
        for iFace = 1:element.numFaces;    
            count = count+1;
            T2D(count,:) = getFace( mesh.T(iElem,:) , element, iFace );
            qualityFaces(count) = mesh.quality(iElem);
            CM(:,count) = xyz;
            if(isfield(mesh,'qualities'))
               for i=1:length(namesQ)
                   qualitiesSubelements.(namesQ{i})(count) = mesh.qualities.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'distortions'))
               for i=1:length(namesQ)
                   distortionsSubelements.(namesQ{i})(count) = mesh.distortions.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'errorDisc'))
               for i=1:length(namesE)
                   errorDisc.(namesE{i})(count) = mesh.errorDisc.(namesE{i})(iElem);
               end
            end
        end
    end
    mesh2D = mesh;
    mesh2D.T = T2D;
    mesh2D.CM = CM;
    mesh2D.quality = qualityFaces;
    if(isfield(mesh,'qualities'))
       mesh2D.qualities = qualitiesSubelements;
    end
    if(isfield(mesh,'distortions'))
       mesh2D.distortions = distortionsSubelements;
    end
    if(isfield(mesh,'errorDisc'))
       for i=1:length(namesE)
           mesh.errorDisc = errorDisc;
       end
    end
    switch(element.type)
        case 'tet'
            element.type = 'tri'; element.dim = 2;
            [ elementTri ] = defineElement(element);
            mesh2D.element = elementTri;
        case 'hex'
            element.type = 'quad'; element.dim = 2;
            [ elementQuad ] = defineElement(element);
            mesh2D.element = elementQuad;
    end
end
