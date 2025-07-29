function [nodalDistortion,meshLin] = computeNodalDisortion(options,mesh)

    options.check = 'edges';
    options.tol = 1e-6;
%     options.check = 'det';
%     options.tol = 1e-2;

    if(isfield(options,'pointDistSrf') && options.pointDistSrf)
        [nodalDistortion,meshLin] = computeSurfaceNodalDisortion(options,mesh);
    else
        [nodalDistortion,meshLin] = computeSolidNodalDisortion(options,mesh);
    end

end

function [nodalDistortion,meshLin] = computeSurfaceNodalDisortion(options,mesh)
    expOrder = options.exportOrder;

    nElem = size(mesh.T,1);
	
    equiElemExpOrder = setDefaulElement('3D','tet',expOrder);
    equiElemExpOrder.distribution = 'equispaced';
    numNodesElemExp = equiElemExpOrder.numNod;
    
    numNewNodes = nElem*numNodesElemExp;

    numFields = 2;    
    nodalDistortion = zeros(numNewNodes,numFields);
    
    %% shapeFunctions 
    [reference_coord]=giveReferencePoints(mesh.element);
    evalPoints = giveReferencePoints(equiElemExpOrder);
    [shapeFunctions]=getShapeFunctions(mesh.element,reference_coord,evalPoints);
    %% reestructure vectorial information
    nNodElem = size(mesh.T,2);
    XX = zeros(nElem,1,size(mesh.X,1),nNodElem);
    for iElem =1:nElem
        XX(iElem,1,:,:)=mesh.X(:,mesh.T(iElem,:));
    end
    %% compute Jacobian
    Dphi = computeDphi_HO(XX,shapeFunctions);
    %% evaluate distortion of the jacobian
    WinvIni = mesh.idealInitial;
    WinvGeo = mesh.idealGeometric;
    delta = 0;
    select_LinearDeviationMeasure = 'algebraicShape';
    
    etaIni = computeLocalEta(...
        Dphi,WinvIni,delta,select_LinearDeviationMeasure,[]);
    etaGeo = computeLocalEta(...
        Dphi,WinvGeo,delta,select_LinearDeviationMeasure,[]);

    
    [i] = find(isnan(etaIni));
    etaIni(i) = 1e6;
    [i] = find(isinf(etaIni));
    etaIni(i) = 1e6;
    [i] = find(isnan(etaGeo));
    etaGeo(i) = 1e6;
    [i] = find(isinf(etaGeo));
    etaGeo(i) = 1e6;
    
    %% restructure info
    meshLin.element = setDefaulElement('2D','tri',1);
    elemExpSrf = setDefaulElement('2D','tri',expOrder);
    refLinMesh = giveLinearMesh(1:elemExpSrf.numNod,elemExpSrf);
    TLinRef = refLinMesh.T;
    numSubelemsElem = size(TLinRef,1);
    numLinElem = nElem*numSubelemsElem*4;
    meshLin.T = zeros(numLinElem,3);
    meshLin.X = zeros(3,numNewNodes);
    
    %---
        CM =zeros(3,numLinElem);
        qualitySubelements = zeros(numLinElem,1);
        if(isfield(mesh,'qualities'))
           namesQ = fieldnames(mesh.qualities);
           for i=1:length(namesQ)
               qualitiesSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
           end
        end
        if(isfield(mesh,'distortions'))
           namesQ = fieldnames(mesh.distortions);
           for i=1:length(namesQ)
               distortionsSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
           end
        end
        if(isfield(mesh,'errorDisc'))
           namesE = fieldnames(mesh.errorDisc);
           for i=1:length(namesE)
               errorDisc.(namesE{i}) = zeros(size(qualitySubelements));
           end
        end
    %---
    
%     mesh
%     [nodesToSmooth,objectLayers] = giveFreeNodesNlayers(...
%        mesh.objectElements(:,1),1,mesh.T,size(mesh.X,2),mesh.element);  
   
%     nodalDistortion(:,1) = 1./reshape(etaIni,size(etaIni,1)*size(etaIni,3),1);
%     nodalDistortion(:,2) = 1./reshape(etaGeo,size(etaIni,1)*size(etaIni,3),1);

    elementSrf = setDefaulElement('2D','tri',mesh.element.order);

    numDegIni = 0;
    numDegGeo = 0;
    listDegEl = [];
    countNode = 0;
    countElem = 0;
    countLin = 0;
    countHO = 0;
    for iElem = 1:nElem
        linEl0Face = countElem+1;
        for iFace = 1:mesh.element.numFaces;
            Tface = getFace( mesh.T(iElem,:) , mesh.element, iFace );
            TFaceRef = getFace( 1:equiElemExpOrder.numNod , equiElemExpOrder, iFace );
            vali = etaIni(iElem,:,TFaceRef);
            valg = etaGeo(iElem,:,TFaceRef);
            if(options.activeCostReduction && ...
                    checkLinearFace(elementSrf,mesh.X(:,Tface)))
                countLin = countLin+1;

                nod0 = countNode + 1;
                nod1 = countNode + 3;
                countNode = nod1;
                linEl0 = countElem+1;
                linEl1 = linEl0;            
                countElem = linEl1;

                nodalDistortion(nod0:nod1,1) = mean(1./vali(:));
                nodalDistortion(nod0:nod1,2) = mean(1./valg(:));

                meshLin.T(linEl0:linEl1,:) = nod0:nod1;
                meshLin.X(:,nod0:nod1) = mesh.X(:,Tface(1:3));
            else
                countHO = countHO +1;
                nod0 = 1 + countNode;
                nod1 = countNode + elemExpSrf.numNod;
                countNode = nod1;
                nodalDistortion(nod0:nod1,1) = 1./vali(:);
                nodalDistortion(nod0:nod1,2) = 1./valg(:);

                linEl0 = 1 + countElem;
                linEl1 = countElem + numSubelemsElem;
                countElem = linEl1;
                Telem = nod0:nod1;
                meshLin.T(linEl0:linEl1,:) = Telem(TLinRef);
                meshLin.X(:,nod0:nod1)= mesh.X(:,mesh.T(iElem,:))*shapeFunctions(:,TFaceRef,1);
            end

%         if(isempty(find(objectLayers == iElem)))
%         nodalDistortion(nod0:nod1,1) = 1;
%         nodalDistortion(nod0:nod1,2) = 1;
%         end
        end
        linEl1Face = countElem;
        numNewElements = linEl1Face-linEl0Face+1;
        
        vali = etaIni(iElem,:,:);
        if(max(max(1./vali<0.05))==1)
            listDegEl = [listDegEl iElem];
        end
        numDegIni = numDegIni + max(max(1./vali<0.05));
        numDegGeo = numDegGeo + max(max(1./valg<0.05));
        
        %----
            if(isfield(mesh,'quality') )
                qualitySubelements(linEl0Face:linEl1Face)=mesh.quality(iElem);
            end
            CM(:,linEl0Face:linEl1Face) = mesh.CM(:,iElem)*ones(1,numNewElements);
            if(isfield(mesh,'qualities'))
               for i=1:length(namesQ)
                   qualitiesSubelements.(namesQ{i})(linEl0Face:linEl1Face) = mesh.qualities.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'distortions'))
               for i=1:length(namesQ)
                   distortionsSubelements.(namesQ{i})(linEl0Face:linEl1Face) = mesh.distortions.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'errorDisc'))
               for i=1:length(namesE)
                   errorDisc.(namesE{i})(linEl0Face:linEl1Face) = mesh.errorDisc.(namesE{i})(iElem);
               end
            end
        %----
    end
    
    if(length(listDegEl)<15)
        listDegEl
    end
    numDegIni
    
    nodalDistortion = nodalDistortion(1:countNode,:);
    meshLin.X = meshLin.X(:,1:countNode);
    meshLin.T = meshLin.T(1:countElem,:);
    %----
        if(isfield(mesh,'quality'))
            meshLin.quality = qualitySubelements(1:countElem);
        end
        meshLin.CM = CM(:,1:countElem);
        if(isfield(mesh,'qualities'))
           for i=1:length(namesQ)
               meshLin.qualities.(namesQ{i}) = qualitiesSubelements.(namesQ{i})(1:countElem);
           end
        end
        if(isfield(mesh,'distortions'))
           for i=1:length(namesQ)
               meshLin.distortions.(namesQ{i}) = distortionsSubelements.(namesQ{i})(1:countElem);
           end
        end
        if(isfield(mesh,'errorDisc'))
           for i=1:length(namesE)
               meshLin.errorDisc.(namesE{i}) = errorDisc.(namesE{i})(1:countElem);
           end
        end
    %----
    fprintf('Exported elements: Total -> Lin / HO  ===> %d -> %d / %d \n',nElem,countLin,countHO)
    meshLin    
    writeNodeStatistics(nodalDistortion,mesh.fileName,numDegIni,numDegGeo);
end

function [nodalDistortion,meshLin] = computeSolidNodalDisortion(options,mesh)

    expOrder = options.exportOrder;

    nElem = size(mesh.T,1);
	
    equiElemExpOrder = setDefaulElement('3D','tet',expOrder);
    equiElemExpOrder.distribution = 'equispaced';
    numNodesElemExp = equiElemExpOrder.numNod;
    
    numNewNodes = nElem*numNodesElemExp;

    numFields = 2;    
    nodalDistortion = zeros(numNewNodes,numFields);
    
    %% shapeFunctions 
    [reference_coord]=giveReferencePoints(mesh.element);
    evalPoints = giveReferencePoints(equiElemExpOrder);
    [shapeFunctions]=getShapeFunctions(mesh.element,reference_coord,evalPoints);
    %% reestructure vectorial information
    nNodElem = size(mesh.T,2);
    XX = zeros(nElem,1,size(mesh.X,1),nNodElem);
    for iElem =1:nElem
        XX(iElem,1,:,:)=mesh.X(:,mesh.T(iElem,:));
    end
    %% compute Jacobian
    Dphi = computeDphi_HO(XX,shapeFunctions);
    %% evaluate distortion of the jacobian
    WinvIni = mesh.idealInitial;
    WinvGeo = mesh.idealGeometric;
    delta = 0;
    if(isfield(options,'distortion'))
        select_LinearDeviationMeasure = options.distortion;
    else
        select_LinearDeviationMeasure = 'algebraicShape';
    end
    
    etaIni = computeLocalEta(...
        Dphi,WinvIni,delta,select_LinearDeviationMeasure,[]);
    etaGeo = computeLocalEta(...
        Dphi,WinvGeo,delta,select_LinearDeviationMeasure,[]);

    [i] = find(isnan(etaIni));
    etaIni(i) = 1e6;
    [i] = find(isinf(etaIni));
    etaIni(i) = 1e6;
    [i] = find(isnan(etaGeo));
    etaGeo(i) = 1e6;
    [i] = find(isinf(etaGeo));
    etaGeo(i) = 1e6;

    %% restructure info
    meshLin.element = setDefaulElement('3D','tet',1);
    refLinMesh = giveLinearMesh(1:equiElemExpOrder.numNod,equiElemExpOrder);
    TLinRef = refLinMesh.T;
    numSubelemsElem = size(TLinRef,1);
    numLinElem = nElem*numSubelemsElem;
    meshLin.T = zeros(numLinElem,4);
    meshLin.X = zeros(3,numNewNodes);
    
    
    %---
        CM =zeros(3,numLinElem);
        qualitySubelements = zeros(numLinElem,1);
        if(isfield(mesh,'qualities'))
           namesQ = fieldnames(mesh.qualities);
           for i=1:length(namesQ)
               qualitiesSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
           end
        end
        if(isfield(mesh,'distortions'))
           namesQ = fieldnames(mesh.distortions);
           for i=1:length(namesQ)
               distortionsSubelements.(namesQ{i}) = zeros(size(qualitySubelements));
           end
        end
        if(isfield(mesh,'errorDisc'))
           namesE = fieldnames(mesh.errorDisc);
           for i=1:length(namesE)
               errorDisc.(namesE{i}) = zeros(size(qualitySubelements));
           end
        end
    %---
    
%     mesh
%     [nodesToSmooth,objectLayers] = giveFreeNodesNlayers(...
%        mesh.objectElements(:,1),1,mesh.T,size(mesh.X,2),mesh.element);  
   
%     nodalDistortion(:,1) = 1./reshape(etaIni,size(etaIni,1)*size(etaIni,3),1);
%     nodalDistortion(:,2) = 1./reshape(etaGeo,size(etaIni,1)*size(etaIni,3),1);
    numDegIni = 0;
    numDegGeo = 0;
    listDegEl = [];
    countNode = 0;
    countElem = 0;
    countLin = 0;
    countHO = 0;
    for iElem = 1:nElem
        if(options.activeCostReduction && ...
                checkLinearElement(mesh.element,mesh.X(:,mesh.T(iElem,:)),options))
            countLin = countLin+1;
            
            nod0 = countNode + 1;
            nod1 = countNode + 4;
            countNode = nod1;
            linEl0 = countElem+1;
            linEl1 = linEl0;            
            countElem = linEl1;

            
            vali = etaIni(iElem,:,:);
            valg = etaGeo(iElem,:,:);
            nodalDistortion(nod0:nod1,1) = mean(1./vali(:));
            nodalDistortion(nod0:nod1,2) = mean(1./valg(:));
            
            meshLin.T(linEl0:linEl1,:) = nod0:nod1;
            meshLin.X(:,nod0:nod1) = mesh.X(:,mesh.T(iElem,1:4));
            numNewElements = 1;
        else
            countHO = countHO +1;
            nod0 = 1 + countNode;
            nod1 = countNode + equiElemExpOrder.numNod;
            countNode = nod1;
            vali = etaIni(iElem,:,:);
            valg = etaGeo(iElem,:,:);
            nodalDistortion(nod0:nod1,1) = 1./vali(:);
            nodalDistortion(nod0:nod1,2) = 1./valg(:);

            linEl0 = 1 + countElem;
            linEl1 = countElem + numSubelemsElem;
            countElem = linEl1;
            Telem = nod0:nod1;
            meshLin.T(linEl0:linEl1,:) = Telem(TLinRef);
            meshLin.X(:,nod0:nod1) = mesh.X(:,mesh.T(iElem,:))*shapeFunctions(:,:,1);
            numNewElements = numSubelemsElem;
        end

%         if(isempty(find(objectLayers == iElem)))
%         nodalDistortion(nod0:nod1,1) = 1;
%         nodalDistortion(nod0:nod1,2) = 1;
%         end

        if(max(max(1./vali<0.05))==1)
            listDegEl = [listDegEl iElem];
        end
        numDegIni = numDegIni + max(max(1./vali<0.05));
        numDegGeo = numDegGeo + max(max(1./valg<0.05));
        
        %----
            if(isfield(mesh,'quality') )
                qualitySubelements(linEl0:linEl1)=mesh.quality(iElem);
            end
            if(isfield(mesh,'CM') )
                CM(:,linEl0:linEl1) = mesh.CM(:,iElem)*ones(1,numNewElements);
            end
            if(isfield(mesh,'qualities'))
               for i=1:length(namesQ)
                   qualitiesSubelements.(namesQ{i})(linEl0:linEl1) = mesh.qualities.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'distortions'))
               for i=1:length(namesQ)
                   distortionsSubelements.(namesQ{i})(linEl0:linEl1) = mesh.distortions.(namesQ{i})(iElem);
               end
            end
            if(isfield(mesh,'errorDisc'))
               for i=1:length(namesE)
                   errorDisc.(namesE{i})(linEl0:linEl1) = mesh.errorDisc.(namesE{i})(iElem);
               end
            end
        %----
    end
    
    if(length(listDegEl)<15)
        listDegEl
    end
    numDegIni
    
    nodalDistortion = nodalDistortion(1:countNode,:);
    meshLin.X = meshLin.X(:,1:countNode);
    meshLin.T = meshLin.T(1:countElem,:);
    %----
        if(isfield(mesh,'quality'))
            meshLin.quality = qualitySubelements(1:countElem);
        end
        if(isfield(mesh,'CM') )
            meshLin.CM = CM(:,1:countElem);
        end
        if(isfield(mesh,'qualities'))
           for i=1:length(namesQ)
               meshLin.qualities.(namesQ{i}) = qualitiesSubelements.(namesQ{i})(1:countElem);
           end
        end
        if(isfield(mesh,'distortions'))
           for i=1:length(namesQ)
               meshLin.distortions.(namesQ{i}) = distortionsSubelements.(namesQ{i})(1:countElem);
           end
        end
        if(isfield(mesh,'errorDisc'))
           for i=1:length(namesE)
               meshLin.errorDisc.(namesE{i}) = errorDisc.(namesE{i})(1:countElem);
           end
        end
    %----
    fprintf('Exported elements: Total -> Lin / HO  ===> %d -> %d / %d \n',nElem,countLin,countHO)
    meshLin    
    writeNodeStatistics(nodalDistortion,mesh.fileName,numDegIni,numDegGeo);
end

