function mesh = computePeriodicMaps(mesh)
% mesh.X: Coordinates (n x 2 matrix for 2D mesh)
% mesh.T: Connectivity (m x p matrix for mesh elements)

% Extract mesh dimensions
xCoords = mesh.X(:, 1);
yCoords = mesh.X(:, 2);

% Define tolerance for node matching
tolerance = 1e-6;

% Find boundary nodes based on the coordinates
% Left boundary: x = min(xCoords)
leftBoundary = find(abs(xCoords - min(xCoords)) < tolerance);

% Right boundary: x = max(xCoords)
rightBoundary = find(abs(xCoords - max(xCoords)) < tolerance);

% Bottom boundary: y = min(yCoords)
bottomBoundary = find(abs(yCoords - min(yCoords)) < tolerance);

% Top boundary: y = max(yCoords)
topBoundary = find(abs(yCoords - max(yCoords)) < tolerance);

% Initialize maps for left-right and top-bottom node correspondences
leftRightMap = containers.Map('KeyType', 'double', 'ValueType', 'double');
topBottomMap = containers.Map('KeyType', 'double', 'ValueType', 'double');

% Match left and right boundary nodes
for i = 1:length(leftBoundary)
    leftNode = leftBoundary(i);
    leftCoord = mesh.X(leftNode, :);

    % Find corresponding right node
    rightNode = rightBoundary(find(abs(mesh.X(rightBoundary, 2) - leftCoord(2)) < tolerance, 1));
    if ~isempty(rightNode)
        leftRightMap(leftNode) = rightNode;
    else
        error('Not found periodic image of the node')
    end
end

% Match bottom and top boundary nodes
for i = 1:length(bottomBoundary)
    bottomNode = bottomBoundary(i);
    bottomCoord = mesh.X(bottomNode, :);

    % Find corresponding top node
    topNode = topBoundary(find(abs(mesh.X(topBoundary, 1) - bottomCoord(1)) < tolerance, 1));
    if ~isempty(topNode)
        topBottomMap(topNode) = bottomNode;
    else
        error('Not found periodic image of the node')
    end
end

% Output the maps
periodic_maps.left  = leftBoundary;
periodic_maps.right = rightBoundary;
periodic_maps.top = topBoundary;
periodic_maps.bottom = bottomBoundary;
periodic_maps.leftRight = leftRightMap;
periodic_maps.topBottom = topBottomMap;
mesh.periodic_maps = periodic_maps;


