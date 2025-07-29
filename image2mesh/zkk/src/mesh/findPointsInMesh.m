function [elements,shapeF,typeCoord,shapeF_byElement]=findPointsInMesh(mesh,points,options)

% first we should use a bin/octree structure to make the routine more
% efficient

%location.element: element where it is
%location.coord: baricentric/parametric coordinates in the element
%location.type: barycentric/parametric...?

jacobianRef = 2;

%global element;
%element = mesh.element;

if(nargin>2 && isfield(options,'classifyByElement'))
    classifyByElement = options.classifyByElement;    
else
    classifyByElement = false;
    classifyByElement = true;
end
if(nargin>2 && isfield(options,'stopWithNotFound'))
    stopWithNotFound = options.stopWithNotFound;    
else
    stopWithNotFound = true;
end
if(nargin>2 && isfield(options,'bin'))
    useBins = options.bin;    
else
    useBins = false;
end

numPoints = size(points,1);
numElements = size(mesh.T,1);

elements = zeros(numPoints,1);
shapeF = zeros(numPoints,3);


if(classifyByElement)
   %maxPointsPerElement = 200000;
   %SHelem = zeros(numElements,maxPointsPerElement,3);
   %pointsElem = zeros(numElements,maxPointsPerElement);
   countPointsInElement = zeros(numElements,1);
   SHelem = cell(numElements,1);
   pointsElem = cell(numElements,1);
else
    shapeF_byElement = [];
end

if(useBins)
    %disp('splitting bins')
    numBinsDirection = options.numBins;
    [bins] = binMesh(mesh,numBinsDirection);
    %disp('done splitting bins')
else
    %?
    numBinsDirection = 0;%???
    error('Now use with bins and change later')
end


%bins.elements = binElements;
%bins.field = binField;
numBins = length(bins.elements);
[loc] = locate(points, bins.field);
fprintf('\n Num bins %d. Currently: ',numBins)
for ibin=1:numBins
    
    fprintf('%d ',ibin)
    
    pointsInBin = find(loc.l==ibin);
    meshInBin = mesh;
    meshInBin.T = mesh.T(bins.elements{ibin},:);
    %meshInBin.T = mesh.T(bins.elements{ibin},[2 3 1 ]);
    
%     figure(1)
%     hold on
%     plot(points(pointsInBin,1),points(pointsInBin,2),'*')
%     pause()
    
    if(~isempty(bins.elements{ibin}))
        [elements(pointsInBin),shapeF(pointsInBin,:),typeCoord,shapeF_byElement_bin]=...
            findPointsInMesh_vec(meshInBin,points(pointsInBin,:),classifyByElement,stopWithNotFound);

        for ielem=1:length(bins.elements{ibin})
            theElem = bins.elements{ibin}(ielem);
            
            newNumPts = shapeF_byElement_bin.numPointsInElem(ielem);
            %newIndex = (countPointsInElement(theElem)+1):(countPointsInElement(theElem)+newNumPts);

            %pointsElem(theElem,newIndex)=pointsInBin(shapeF_byElement_bin.points(ielem,1:newNumPts));
            pointsElem{theElem}=[pointsElem{theElem}; pointsInBin(shapeF_byElement_bin.points{ielem})]; 
            
            SHelem{theElem} = [ SHelem{theElem}; shapeF_byElement_bin.shapeF{ielem}*jacobianRef];
            
            countPointsInElement(theElem) = countPointsInElement(theElem)+newNumPts;
        end
        
    end
end

fprintf('Num of found elements: %d (out of %d)\n',sum(countPointsInElement),numPoints)

% if(classifyByElement)
%    maxPointsPerElement = 1000;
%    SHelem = zeros(numElements,maxPointsPerElement,3);
%    pointsElem = zeros(numElements,maxPointsPerElement);
%    countPointsInElement = zeros(numElements,1);
%    
%     totalMaxPointsInElem = 0;
%     for ielem = 1:numElements
%        pointsInElem = find(elements==ielem); 
%        numPointsInElem = length(pointsInElem);
%        totalMaxPointsInElem = max([totalMaxPointsInElem numPointsInElem]);
%        
%        SHelem(ielem,1:numPointsInElem,:) = shapeF(pointsInElem,:);
%        pointsElem(ielem,1:numPointsInElem) = pointsInElem;
%        countPointsInElement(ielem) = numPointsInElem;
%     end
%     SHelem = SHelem(:,1:totalMaxPointsInElem,:);
%     pointsElem = pointsElem(:,1:totalMaxPointsInElem);
% end
% 
if(classifyByElement)
    shapeF_byElement.numPointsInElem = countPointsInElement;
    shapeF_byElement.shapeF = SHelem;%(:,1:max(countPointsInElement),:);
    shapeF_byElement.points = pointsElem;
end 

end

function [elements,shapeF,typeCoord,shapeF_byElement]=findPointsInMesh_vec(mesh,points,classifyByElement,stopWithNotFound)
shapeF_byElement = [];

typeCoord = 'barycentric';

numPoints = size(points,1);
numElements = size(mesh.T,1);

if(classifyByElement)
   maxPointsPerElement = 1000;
   SHelem = cell(numElements,1);%zeros(numElements,maxPointsPerElement,3);
   pointsElem = cell(numElements,1);%zeros(numElements,maxPointsPerElement);
   countPointsInElement = zeros(numElements,1);
else
    shapeF_byElement = [];
end

xcm = (mesh.X(mesh.T(:,1),1)+mesh.X(mesh.T(:,2),1)+mesh.X(mesh.T(:,3),1))/3.0;
ycm = (mesh.X(mesh.T(:,1),2)+mesh.X(mesh.T(:,2),2)+mesh.X(mesh.T(:,3),2))/3.0;

CM = zeros(numElements,1,2);
CM(:,:,1) = xcm;
CM(:,:,2) = ycm;

points_vec = reshape(points,[1 size(points)]);

maxNeighSearch = size(mesh.T,1);%50;

elements = zeros(numPoints,1);
shapeF = zeros(numPoints,3);

dist = bsxfun(@minus,CM(:,1,:),points_vec(1,:,:));    
dist = sum(dist.*dist,3);

if(numElements>1)
    [mind,closestElem] = min(dist);
else
    closestElem = ones(numPoints,1);
end
Xelems = zeros(length(closestElem),3,2);
Xelems(:,:,1) = [mesh.X(mesh.T(closestElem,1),1) mesh.X(mesh.T(closestElem,2),1) mesh.X(mesh.T(closestElem,3),1)];
Xelems(:,:,2) = [mesh.X(mesh.T(closestElem,1),2) mesh.X(mesh.T(closestElem,2),2) mesh.X(mesh.T(closestElem,3),2)];
[isInElement,shcomputed]=findPointInTriangle_vec(points,Xelems);
remainingPoints = find(~isInElement);
foundPoints = find(isInElement);
elements(foundPoints) = closestElem(foundPoints);
shapeF(foundPoints,:) = shcomputed(foundPoints,:);

closestElem = closestElem(remainingPoints);

numSearch = 0;
while(~isempty(remainingPoints) && numSearch < maxNeighSearch && numElements>1)
    
    linInd = sub2ind(size(dist),closestElem',remainingPoints);
    dist(linInd) = NaN;
    
    [mind,closestElem] = min(dist(:,remainingPoints));
    Xelems = zeros(length(closestElem),3,2);
    Xelems(:,:,1) = [mesh.X(mesh.T(closestElem,1),1) mesh.X(mesh.T(closestElem,2),1) mesh.X(mesh.T(closestElem,3),1)];
    Xelems(:,:,2) = [mesh.X(mesh.T(closestElem,1),2) mesh.X(mesh.T(closestElem,2),2) mesh.X(mesh.T(closestElem,3),2)];
    [isInElement,shcomputed]=findPointInTriangle_vec(points(remainingPoints,:),Xelems);
    locrem = find(~isInElement);
    locfound = find(isInElement);

    foundPoints = remainingPoints(locfound);
    elements(foundPoints) = closestElem(locfound);
    shapeF(foundPoints,:) = shcomputed(locfound,:);
    
    closestElem = closestElem(locrem);
    remainingPoints = remainingPoints(locrem);

    numSearch = numSearch+1;
end

if(stopWithNotFound)
%     if(max(find(theElement==0)))%~isempty(remainingPoints) && numSearch == maxNeighSearch ) %isempty(dist))
    if(max(find(closestElem==0)))%~isempty(remainingPoints) && numSearch == maxNeighSearch ) %isempty(dist))
        disp(remainingPoints);
        disp(points(remainingPoints,:));
        error('Not found point')
    else
        %theElement(remainingPoints) = 0;
        %sh(remainingPoints,:) = 0;
    end
end

if(classifyByElement)
   %maxPointsPerElement = 1000;
   %SHelem = zeros(numElements,maxPointsPerElement,3);
   %pointsElem = zeros(numElements,maxPointsPerElement);
   %countPointsInElement = zeros(numElements,1);
   
    totalMaxPointsInElem = 0;
    for ielem = 1:numElements
       pointsInElem = find(elements==ielem); 
       numPointsInElem = length(pointsInElem);
       totalMaxPointsInElem = max([totalMaxPointsInElem numPointsInElem]);
       
       %SHelem(ielem,1:numPointsInElem,:) = shapeF(pointsInElem,:);
       %pointsElem(ielem,1:numPointsInElem) = pointsInElem;
       countPointsInElement(ielem) = numPointsInElem;
       SHelem{ielem} =  shapeF(pointsInElem,:);
       pointsElem{ielem} = pointsInElem;
    end
    %SHelem = SHelem(:,1:totalMaxPointsInElem,:);
    %pointsElem = pointsElem(:,1:totalMaxPointsInElem);
end

if(classifyByElement)
    shapeF_byElement.numPointsInElem = countPointsInElement;
    shapeF_byElement.shapeF = SHelem;%(:,1:max(countPointsInElement),:);
    shapeF_byElement.points = pointsElem;
end 

end


function [isInElement,shapeF]=findPointInTriangle_vec(points,Xelems)

%    global element;

    tolZero = 1e-12;

    % Equivalent ways
    x = points(:,1);
    y = points(:,2);
    x1 = Xelems(:,1,1);x2 = Xelems(:,2,1);x3 = Xelems(:,3,1);
    y1 = Xelems(:,1,2);y2 = Xelems(:,2,2);y3 = Xelems(:,3,2);
    lambda1    = ( (y2-y3).*(x-x3)  + (x3-x2).*(y-y3) )./...
                 ( (y2-y3).*(x1-x3) + (x3-x2).*(y1-y3) );
    lambda2    = ( (y3-y1).*(x-x3)  + (x1-x3).*(y-y3) )./...
                 ( (y2-y3).*(x1-x3) + (x3-x2).*(y1-y3) );
    lambda3    = 1.0 - lambda1 - lambda2;
    shapeF = [lambda1 lambda2 lambda3];
            
    % Equivalent ways
%     A = ones(size(Xelems,1),3,3);
%     A(:,1:2,1) = Xelems(:,1,:);
%     A(:,1:2,2) = Xelems(:,2,:);
%     A(:,1:2,3) = Xelems(:,3,:);
%     detA = determinantNDim(A);
%     Ainv = invertMultidimensionalMatrix(A,detA);
%     lambda1 = Ainv(:,1,1).*points(:,1) + Ainv(:,1,2).*points(:,2) + Ainv(:,1,3);
%     lambda2 = Ainv(:,2,1).*points(:,1) + Ainv(:,2,2).*points(:,2) + Ainv(:,2,3);
%     lambda3 = Ainv(:,3,1).*points(:,1) + Ainv(:,3,2).*points(:,2) + Ainv(:,3,3);
%     shapeF = [lambda1 lambda2 lambda3];
    
    %isInElement = min(shapeF')>=-tolZero ;%&& max(coords)<=1;
    %isInElement = and(and(lambda1(:)>0,(lambda2(:)>0)), (lambda3(:)>0));
    
    isInElement = min(shapeF,[],2)>=-tolZero ;%&& max(coords)<=1;

    %% compute the desired shape functions
%     elemCoord = [-1 -1; 1  -1; -1 1];
%                   
%     chieta = [shapeF*elemCoord(:,1) shapeF*elemCoord(:,2) ];
%     [shapeFunctions]=getShapeFunctions(element,elemCoord,chieta);
%     
%     shapeF(:,:) = shapeFunctions(:,:,1)';
%     
%error('map to the reference element and use the shape funcitons!!!!')
    
end

function [bins] = binMesh(mesh,numBinsDirection)

    x = mesh.X(:,1);
    y = mesh.X(:,2);
    
    myEps = 1e-10;
    minx = min(x) - myEps; 
    maxx = max(x) + myEps;
    miny = min(y) - myEps; 
    maxy = max(y) + myEps;

    lx = maxx-minx;
    ly = maxy-miny;
    
    if(lx>=ly)
        nx = numBinsDirection;
        hx = lx/nx;
        ny = ceil(ly/hx);
        hy = ly/ny;
    else
        ny = numBinsDirection;
        hy = ly/ny;
        nx = ceil(lx/hy);
        hx = lx/nx;
    end
    
    vertices(:,1,:) = mesh.X(mesh.T(:, 1 ),:);
    vertices(:,2,:) = mesh.X(mesh.T(:, 2 ),:);
    vertices(:,3,:) = mesh.X(mesh.T(:, 3 ),:);
    vertices(:,4,:) = mesh.X(mesh.T(:, 1 ),:);
    
    binField.hx = hx;
    binField.hy = hy;
    binField.m = nx;%+1;
    binField.n = ny;%+1;
    binField.x0 = [minx miny];
    
    numBins = binField.m*binField.n;%(nx+1)*(ny+1) ;%nx*ny;
    
%     xbin = minx:hx:(minx+ hx*binField.m);
%     ybin = miny:hy:(miny+ hy*binField.n);
%     [xx, yy] = meshgrid(xbin,ybin);
%     figure(1)
%     hold on
%     plot(xx(:),yy(:),'ko')
    
    %binField.m = nx+1;
    %binField.n = ny+1;
    
%     firstNode  = mesh.X(mesh.T(:,1),:);
%     secondNode = mesh.X(mesh.T(:,2),:);
%     thirdNode  = mesh.X(mesh.T(:,3),:);
%     [loc1] = locate(firstNode, binField);
%     [loc2] = locate(secondNode, binField);
%     [loc3] = locate(thirdNode, binField);
%     
%     binElements = cell(numBins,1);
%     for ibin = 1:numBins
%        binElements{ibin} = unique([find(loc1.l==ibin)  ; find(loc2.l==ibin) ; find(loc3.l==ibin)]);
%     end

    binElements = cell(numBins,1);
    
    numSamplings = 100;
    %n1 = zeros(size(mesh.T,1),size(mesh.X,2));
    %n2 = zeros(size(mesh.T,1),size(mesh.X,2));
    for iedge=1:size(mesh.T,2)
        n1(:,1) = vertices(:,iedge,1);
        n1(:,2) = vertices(:,iedge,2);
        n2(:,1) = vertices(:,iedge+1,1);
        n2(:,2) = vertices(:,iedge+1,2);
        
        h = 1.0/numSamplings;
        v = n2-n1;
        for isample=1:numSamplings
            samplePoint = n1 + v*(h*(isample-1));
            [loc] = locate(samplePoint, binField);
            for ibin = 1:numBins
                binElements{ibin} = unique([binElements{ibin}; find(loc.l==ibin)]);
            end
        end
    end
    
    binInElement = zeros(numBins,1);
    for ibin = 1:numBins
        if(isempty(binElements{ibin}))
            %disp('Empty bin')
            [i, j] = ind2sub([binField.m binField.n ],ibin);
            hx = binField.hx;
            hy = binField.hy;
            x0 = binField.x0;    
            thePoint = x0 + [(i-1)*hx , (j-1)*hy];
            
            elemContainingBin=findPointsInMesh_vec(mesh,thePoint,false,false);
            if(elemContainingBin>0)
                binElements{ibin} = elemContainingBin;
                binInElement(ibin) = 1;
            end
            %elemContainingBin
        end
    end
    
    bins.elements = binElements;
    bins.field = binField;
    bins.inElement = binInElement;
    
%     numDivisions = 10;
%     hedge = 1.0/numDivisions;
%     chi = 0:hedge:1;
%     %buscar tots els numDivisions punts a cada edge
%     for iedge = 1:3
%        
%         
%     end
end



%%



function [isInElement,shapeF]=findPointInTriangle(point,X)

    tolZero = 1e-14;

%     v1 = X(2,:)-X(1,:);
%     v2 = X(3,:)-X(2,:);
%     v3 = X(1,:)-X(3,:);
%     n1 = [-v1(2); v1(1)];
%     n2 = [-v2(2); v2(1)];
%     n3 = [-v3(2); v3(1)];
%     sign1 = (point-X(1,:))*n1 ;
%     sign2 = (point-X(2,:))*n2 ;
%     sign3 = (point-X(3,:))*n3 ;
%     
%     isInElement = sign1>= -tolZero && sign2>= -tolZero && sign3>= -tolZero;
%     isInElementOld = isInElement;

%     A = [ X(1,:)-X(2,:) ; X(1,:)-X(3,:)]';
%     lambda23 = A\(X(1,:)-point)';
%     lambda1 = 1 - sum(lambda23);
%     
%     shapeF = [ lambda1 lambda23(1) lambda23(2) ];
%     isInElement = (lambda23(1)>=0 && lambda23(1)<=1 && lambda23(2)>=0 && lambda23(2)<=1 && lambda1>=0 && lambda1<=1);

    x = point(1);
    y = point(2);
    x1 = X(1,1);x2 = X(2,1);x3 = X(3,1);
    y1 = X(1,2);y2 = X(2,2);y3 = X(3,2);
    lambda1 = ( (y2-y3)*(x-x3)  + (x3-x2)*(y-y3) )/...
              ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) );
    lambda2 = ( (y3-y1)*(x-x3)  + (x1-x3)*(y-y3) )/...
              ( (y2-y3)*(x1-x3) + (x3-x2)*(y1-y3) );
    lambda3 = 1 - lambda1 - lambda2;
    
    shapeF = [lambda1 lambda2 lambda3];
          
    isInElement = min(shapeF)>=-tolZero ;%&& max(coords)<=1;

%     if(isInElement ~= isInElementOld)
%         disp('m....')
%     end
end




function [elements,shapeF,typeCoord,shapeF_byElement]=findPointsInMesh_bins_sec(mesh,points,options)

% first we should use a bin/octree structure to make the routine more
% efficient

%location.element: element where it is
%location.coord: baricentric/parametric coordinates in the element
%location.type: barycentric/parametric...?

if(nargin>2 && isfield(options,'classifyByElement'))
    classifyByElement = options.classifyByElement;    
else
    classifyByElement = false;
end
if(nargin>2 && isfield(options,'stopWithNotFound'))
    stopWithNotFound = options.stopWithNotFound;    
else
    stopWithNotFound = true;
end
if(nargin>2 && isfield(options,'bin'))
    useBins = options.bin;    
else
    useBins = false;
end

if(useBins)
    disp('splitting bins')
    numBinsDirection = 300;%500;%options.numBins;
    [bins] = binMesh(mesh,numBinsDirection);
    disp('done splitting bins')
end

numPoints = size(points,1);
elements = zeros(numPoints,1);
shapeF = zeros(numPoints,3);
typeCoord = 'barycentric';

if(classifyByElement)
   maxPointsPerElement = 1000;
   SHelem = zeros(size(mesh.T,1),maxPointsPerElement,3);
   pointsElem = zeros(size(mesh.T,1),maxPointsPerElement);
   countPointsInElement = zeros(size(mesh.T,1),1);
else
    shapeF_byElement = [];
end

numElements = size(mesh.T,1);

xcm = (mesh.X(mesh.T(:,1),1)+mesh.X(mesh.T(:,2),1)+mesh.X(mesh.T(:,3),1))/3.0;
ycm = (mesh.X(mesh.T(:,1),2)+mesh.X(mesh.T(:,2),2)+mesh.X(mesh.T(:,3),2))/3.0;

CM = [xcm ycm];

maxNeighSearch = 10;

theElement = 0;

for ipoint = 1:numPoints
    
    if(mod(ipoint,1000)==0)
        disp([ipoint numPoints])
    end
    
    if(useBins)
        [loc] = locate(points(ipoint,:), bins.field);
        
        if(loc.i<1 || loc.i> bins.field.m || loc.j<1 || loc.j> bins.field.n )
            binElements = []; 
            theElement = 0;
            sh = 0;
        else
            ibin = loc.l;
            if(ibin>length(bins.elements))
                loc.i
                loc.j
                loc.l
                bins.field
            end
            binElements = bins.elements{ibin};
        end
    else
        binElements = 1:numElements;
    end
    
    isInElement = false;     
    
    dist = bsxfun(@minus,CM(binElements,:),points(ipoint,:));    
    dist = sum(dist.*dist,2);

    listElements = binElements;
    elem = [];
    numSearch = 0;
    while(~isInElement && ~isempty(dist(1:(end-1))) && numSearch < maxNeighSearch )
        dist(elem) = [];
        listElements(elem) = [];

        [mind,elem] = min(dist);
        theElement = listElements(elem);
        [isInElement,sh]=findPointInTriangle(points(ipoint,:),mesh.X(mesh.T(theElement,:),1:2));

        numSearch = numSearch+1;
    end
    
    if(isempty(dist(1:(end-1))) && numSearch == maxNeighSearch ) %isempty(dist))
        if(stopWithNotFound)
            disp(ipoint);
            disp(points(ipoint,:));
            error('Not found point')
        else
            theElement = 0;
            sh = 0;
        end
    end
    
    % change coords to [-1,1][-1,1][-1,1] triangle ENSURE THIS IS RIGHT
    %coords = coords /2.0;
    
    elements(ipoint) = theElement;
    shapeF(ipoint,:) = sh;
    
    if(classifyByElement)
        if(theElement > 0)
            countPointsInElement(theElement) = countPointsInElement(theElement)+1;
            SHelem(theElement,countPointsInElement(theElement),:) = sh;
            pointsElem(theElement,countPointsInElement(theElement)) = ipoint;
        end
    end
end

if(classifyByElement)
    shapeF_byElement.numPointsInElem = countPointsInElement;
    shapeF_byElement.shapeF = SHelem(:,1:max(countPointsInElement),:);
    shapeF_byElement.points = pointsElem;
end 

end
