function [elements,shapeF,typeCoord,shapeF_byElement]=findPointsInMesh(mesh,points,options)

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

numPoints = size(points,1);
numElements = size(mesh.T,1);

typeCoord = 'barycentric';

if(classifyByElement)
   maxPointsPerElement = 1000;
   SHelem = zeros(numElements,maxPointsPerElement,3);
   pointsElem = zeros(numElements,maxPointsPerElement);
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

maxNeighSearch = 5;

elements = zeros(numPoints,1);
shapeF = zeros(numPoints,3);

dist = bsxfun(@minus,CM(:,1,:),points_vec(1,:,:));    
dist = sum(dist.*dist,3);

[mind,closestElem] = min(dist);
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
while(~isempty(remainingPoints) && numSearch < maxNeighSearch )
    numSearch
    %dist(elem(remainingPoints),remainingPoints) = NaN;
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
    if(max(find(theElement==0)))%~isempty(remainingPoints) && numSearch == maxNeighSearch ) %isempty(dist))
        disp(remainingPoints);
        disp(points(remainingPoints,:));
        error('Not found point')
    else
        %theElement(remainingPoints) = 0;
        %sh(remainingPoints,:) = 0;
    end
end


if(classifyByElement)
   maxPointsPerElement = 1000;
   SHelem = zeros(numElements,maxPointsPerElement,3);
   pointsElem = zeros(numElements,maxPointsPerElement);
   countPointsInElement = zeros(numElements,1);
   
    totalMaxPointsInElem = 0;
    for ielem = 1:numElements
       pointsInElem = find(elements==ielem); 
       numPointsInElem = length(pointsInElem);
       totalMaxPointsInElem = max([totalMaxPointsInElem numPointsInElem]);
       
       SHelem(ielem,1:numPointsInElem,:) = shapeF(pointsInElem,:);
       pointsElem(ielem,1:numPointsInElem) = pointsInElem;
       countPointsInElement(ielem) = numPointsInElem;
    end
    SHelem = SHelem(:,1:totalMaxPointsInElem,:);
    pointsElem = pointsElem(:,1:totalMaxPointsInElem);
end

if(classifyByElement)
    shapeF_byElement.numPointsInElem = countPointsInElement;
    shapeF_byElement.shapeF = SHelem;%(:,1:max(countPointsInElement),:);
    shapeF_byElement.points = pointsElem;
end 

end

function [elements,shapeF,typeCoord,shapeF_byElement]=findPointsInMesh_bins(mesh,points,options)

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

function [isInElement,shapeF]=findPointInTriangle(point,X)

    tolZero = 1e-12;

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

function [isInElement,shapeF]=findPointInTriangle_vec(points,Xelems)

    tolZero = 1e-12;

    x = points(:,1);
    y = points(:,2);
    x1 = Xelems(:,1,1);x2 = Xelems(:,2,1);x3 = Xelems(:,3,1);
    y1 = Xelems(:,1,2);y2 = Xelems(:,2,2);y3 = Xelems(:,3,2);
    lambda1    = ( (y2-y3).*(x-x3)  + (x3-x2).*(y-y3) )./...
                 ( (y2-y3).*(x1-x3) + (x3-x2).*(y1-y3) );
    lambda2    = ( (y3-y1).*(x-x3)  + (x1-x3).*(y-y3) )./...
                 ( (y2-y3).*(x1-x3) + (x3-x2).*(y1-y3) );
    lambda3    = 1 - lambda1 - lambda2;
    
    shapeF = [lambda1 lambda2 lambda3];
          
    %isInElement = min(shapeF')>=-tolZero ;%&& max(coords)<=1;
    isInElement = min(shapeF,[],2)>=-tolZero ;%&& max(coords)<=1;

end

function [bins] = binMesh(mesh,numBinsDirection)

    x = mesh.X(:,1);
    y = mesh.X(:,2);
    
    myEps = 1;
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

    numBins = (nx+1)*(ny+1) ;%nx*ny;
    
%     corners = zeros(nx+1,ny+1);
%      [xcorners, ycorners] = meshgrid(minx:hx:maxx,miny:hy:maxy);
%      size(xcorners')
%     xcorners = xcorners';
%     ycorners = ycorners';
    
    %vertices = mesh.X(mesh.T(:,[1 2 3 1 ]),:);
    firstNode  = mesh.X(mesh.T(:,1),:);
    secondNode = mesh.X(mesh.T(:,2),:);
    thirdNode  = mesh.X(mesh.T(:,3),:);
    
    binField.hx = hx;
    binField.hy = hy;
    binField.m = nx+1;
    binField.n = ny+1;
    binField.x0 = [minx miny];
    
    [loc1] = locate(firstNode, binField);
    [loc2] = locate(secondNode, binField);
    [loc3] = locate(thirdNode, binField);
    
    binElements = cell(numBins,1);
    for ibin = 1:numBins
       binElements{ibin} = unique([find(loc1.l==ibin)  ; find(loc2.l==ibin) ; find(loc3.l==ibin)]);
    end
    
    bins.elements = binElements;
    bins.field = binField;
%     numDivisions = 10;
%     hedge = 1.0/numDivisions;
%     chi = 0:hedge:1;
%     %buscar tots els numDivisions punts a cada edge
%     for iedge = 1:3
%        
%         
%     end
end





