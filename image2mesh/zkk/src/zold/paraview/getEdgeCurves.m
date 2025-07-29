function [mesh1D] = getEdgeCurves(mesh,options)


    if(isfield(options,'exportOrder'))
        factorIncreaseOrder = options.exportOrder / mesh.element.order;
    else
        factorIncreaseOrder = options.factorIncreaseOrder;
    end

    activeCostReduction = options.activeCostReduction;
    
    X=mesh.X; T=mesh.T; element=mesh.element; CM = mesh.CM;
    if(isfield(mesh,'quality'))
        Q = mesh.quality;
    else
        Q = ones(size(mesh.T,1),1);
    end
    [p_order]=element.order;
    feketeDistribution = strcmp(element.distribution,'fekete');
    [shapeFunctions] = getHOcurve_preComputations(...
        p_order,0,feketeDistribution,factorIncreaseOrder);
    newNumPointsEdge=size(shapeFunctions,2);
    X1D = zeros(3,newNumPointsEdge*element.numEdges*size(T,1));
    CM1D = zeros(3,newNumPointsEdge*element.numEdges*size(T,1));
    T1D = zeros(newNumPointsEdge*element.numEdges*size(T,1),2);
    Q1D = zeros(newNumPointsEdge*element.numEdges*size(T,1),1);
    if(isfield(mesh,'qualities'))
       namesQ = fieldnames(mesh.qualities);
       for i=1:length(namesQ)
           qualitiesSubelements.(namesQ{i}) = zeros(size(Q1D));
       end
    end
    if(isfield(mesh,'distortions'))
       namesQ = fieldnames(mesh.distortions);
       for i=1:length(namesQ)
           distortionsSubelements.(namesQ{i}) = zeros(size(Q1D));
       end
    end
    if(isfield(mesh,'errorDisc'))
       namesE = fieldnames(mesh.errorDisc);
       for i=1:length(namesE)
           errorDisc.(namesE{i}) = zeros(size(Q1D));
       end
    end
    countX = 0;
    countT = 0;
    for i=1:size(T,1)
        for iEdge = 1:element.numEdges
            countT = countT+1; countX = countX+1;
            Tedge = getEdge(T(i,:),element,iEdge);
            Xedge = X(:,Tedge);
            if(activeCostReduction && checkLinearEdge(Xedge))
%             if( checkLinearEdge(Xedge) )
                physicalPoints = Xedge(:,[1 length(Tedge)]);
%             countX    
%             countX-1+ size(physicalPoints,2)
%             
%             countT
%              countT-1+ size(physicalPoints,2)-1
            else
                physicalPoints = Xedge*shapeFunctions;
            end
            countXn = countX-1+ size(physicalPoints,2);
            X1D(:,countX:countXn)= physicalPoints;
            countTn = countT-1+ size(physicalPoints,2)-1;
            T1D(countT:countTn,:)= [ (countX:(countXn-1))' ((countX+1):countXn)'];
            CM1D(:,countT:countTn) = CM(:,i)*ones(1,length(countT:countTn));
            Q1D(countT:countTn)= Q(i);
           if(isfield(mesh,'qualities'))
               for ii=1:length(namesQ)
                   qualitiesSubelements.(namesQ{ii})(countT:countTn) = mesh.qualities.(namesQ{ii})(i);
               end
           end
           if(isfield(mesh,'distortions'))
               for ii=1:length(namesQ)
                   distortionsSubelements.(namesQ{ii})(countT:countTn) = mesh.distortions.(namesQ{ii})(i);
               end
           end
           if(isfield(mesh,'errorDisc'))
               for ii=1:length(namesE)
                   errorDisc.(namesE{ii})(countT:countTn) = mesh.errorDisc.(namesE{ii})(i);
               end
           end
            countT = countTn;
            countX = countXn;
        end
    end
    X1D = X1D(:,1:countXn);
    T1D = T1D(1:countTn,:);
    CM1D = CM1D(:,1:countTn);
    Q1D = Q1D(1:countTn);
    if(isfield(mesh,'qualities'))
       for ii=1:length(namesQ)
           qualitiesSubelements.(namesQ{ii}) = qualitiesSubelements.(namesQ{ii})(1:countTn);
       end
    end
    if(isfield(mesh,'distortions'))
       for ii=1:length(namesQ)
           distortionsSubelements.(namesQ{ii}) = distortionsSubelements.(namesQ{ii})(1:countTn);
       end
    end
    if(isfield(mesh,'errorDisc'))
       for ii=1:length(namesE)
           errorDisc.(namesE{ii}) = errorDisc.(namesE{ii})(1:countTn);
       end
    end 
    mesh1D.X = X1D;
    mesh1D.T = T1D;
    mesh1D.CM = CM1D;
    mesh1D.quality = Q1D;
    if(isfield(mesh,'qualities'))
       mesh1D.qualities = qualitiesSubelements;
    end
    if(isfield(mesh,'distortions'))
       mesh1D.distortions = distortionsSubelements;
    end
    if(isfield(mesh,'errorDisc'))
       mesh1D.errorDisc = errorDisc;
    end
    
end



%    if(isfield(mesh,'qualities'))
%        qualitiesSubelements.IdGeo =  zeros(size(Q1D));
%        qualitiesSubelements.IdIni =  zeros(size(Q1D));
%        qualitiesSubelements.ScJac =  zeros(size(Q1D));
%     end
% 
%         if(isfield(mesh,'qualities'))
%            qualitiesSubelements.IdGeo(countT:countTn) =  mesh.qualities.IdGeo(i);
%            qualitiesSubelements.IdIni(countT:countTn) =  mesh.qualities.IdIni(i);
%            qualitiesSubelements.ScJac(countT:countTn) =  mesh.qualities.ScJac(i);
%         end
% 
%    if(isfield(mesh,'qualities'))
%        qualitiesSubelements.IdGeo =  qualitiesSubelements.IdGeo(1:countTn);
%        qualitiesSubelements.IdIni =  qualitiesSubelements.IdIni(1:countTn);
%        qualitiesSubelements.ScJac =  qualitiesSubelements.ScJac(1:countTn);
%    end
% 
%     if(isfield(mesh,'qualities'))
%        mesh1D.qualities = qualitiesSubelements;
%     end