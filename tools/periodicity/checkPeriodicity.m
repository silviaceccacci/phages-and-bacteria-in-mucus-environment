function isPeriodic = checkPeriodicity(mesh)
    % mesh.X: Coordinates (n x 2 matrix for 2D mesh)
    % mesh.T: Connectivity (m x p matrix for mesh elements)
    
    % Extract mesh dimensions
    xCoords = mesh.X(:, 1);
    yCoords = mesh.X(:, 2);
    
    % Find boundary nodes based on the coordinates
    % Define tolerance for node matching
    tolerance = 1e-6;
    
    % Left boundary: x = min(xCoords)
    leftBoundary = find(abs(xCoords - min(xCoords)) < tolerance);
    
    % Right boundary: x = max(xCoords)
    rightBoundary = find(abs(xCoords - max(xCoords)) < tolerance);
    
    % Bottom boundary: y = min(yCoords)
    bottomBoundary = find(abs(yCoords - min(yCoords)) < tolerance);
    
    % Top boundary: y = max(yCoords)
    topBoundary = find(abs(yCoords - max(yCoords)) < tolerance);
    
    % Check if left-right boundaries match (periodicity)
    leftNodes = mesh.X(leftBoundary, :);
    rightNodes = mesh.X(rightBoundary, :);
    
    rightNodes = sort(rightNodes(:,2));
    leftNodes  = sort(leftNodes(:,2));

    % Right nodes should match left nodes with some tolerance
    rightMatch = all(abs(leftNodes - rightNodes) < tolerance);
    
    % Check if bottom-top boundaries match (periodicity)
    bottomNodes = mesh.X(bottomBoundary, :);
    topNodes = mesh.X(topBoundary, :);
    
    bottomNodes = sort(bottomNodes(:,1));
    topNodes  = sort(topNodes(:,1));

    % Top nodes should match bottom nodes with some tolerance
    topMatch = all(abs(bottomNodes - topNodes) < tolerance, 2);

    % If both left-right and bottom-top matches are true, the mesh is periodic
    isPeriodic = all(rightMatch) && all(topMatch);

    if isPeriodic
        disp('The mesh is periodic.');
    else
        disp('The mesh is not periodic.');
        error('Not OK without periodicity.')
    end
end
