function periodicMesh = makePeriodicMeshStokes(mesh, hmin)

    X = mesh.X;
    T = mesh.T;  % keep the existing BAMG connectivity
    
    x0 = min(X(:,1));  y0 = min(X(:,2));
    x1 = max(X(:,1));  y1 = max(X(:,2));
    Lx = x1 - x0;     Ly = y1 - y0;

    X(:,1) = X(:,1) - x0;
    X(:,2) = X(:,2) - y0;

    tol  = 1e-10;
    rtol = hmin * 0.5;

    % Identify boundary nodes
    leftNodes   = find(abs(X(:,1))       < tol);
    rightNodes  = find(abs(X(:,1) - Lx)  < tol);
    bottomNodes = find(abs(X(:,2))       < tol);
    topNodes    = find(abs(X(:,2) - Ly)  < tol);

    % Add missing nodes on each side so left/right and top/bottom match
    newLeftY    = setdiff(X(rightNodes,  2), X(leftNodes,   2), 'stable');
    newRightY   = setdiff(X(leftNodes,   2), X(rightNodes,  2), 'stable');
    newBottomX  = setdiff(X(topNodes,    1), X(bottomNodes, 1), 'stable');
    newTopX     = setdiff(X(bottomNodes, 1), X(topNodes,    1), 'stable');

    nOld = size(X, 1);

    newNodes = [zeros(length(newLeftY),  1),  newLeftY;
                Lx*ones(length(newRightY),1),  newRightY;
                newBottomX, zeros(length(newBottomX), 1);
                newTopX,    Ly*ones(length(newTopX),  1)];

    X = [X; newNodes];

    % For each new node, find the nearest old node and remap T
    % (new boundary nodes land exactly on existing edges — snap them)
    if ~isempty(newNodes)
        for k = 1:size(newNodes, 1)
            newIdx = nOld + k;
            dists  = sum((X(1:nOld,:) - X(newIdx,:)).^2, 2);
            [~, nearest] = min(dists);
            % Remap: wherever T points to nearest, also accept newIdx
            % Replace newIdx references in T with nearest (collapse duplicates)
            T(T == nearest) = newIdx;
        end
    end

    % Remove duplicate/too-close nodes on each boundary (same logic as before)
    leftNodes   = find(abs(X(:,1))      < tol);
    rightNodes  = find(abs(X(:,1) - Lx) < tol);
    bottomNodes = find(abs(X(:,2))      < tol);
    topNodes    = find(abs(X(:,2) - Ly) < tol);

    id_toRemove = [];
    [p, idx] = sort(X(leftNodes, 2));
    for i = 2:length(p)
        if p(i) - p(i-1) < rtol
            id_toRemove(end+1) = leftNodes(idx(i-1));
        end
    end
    [p, idx] = sort(X(rightNodes, 2));
    for i = 2:length(p)
        if p(i) - p(i-1) < rtol
            id_toRemove(end+1) = rightNodes(idx(i-1));
        end
    end
    [p, idx] = sort(X(bottomNodes, 1));
    for i = 2:length(p)
        if p(i) - p(i-1) < rtol
            id_toRemove(end+1) = bottomNodes(idx(i-1));
        end
    end
    [p, idx] = sort(X(topNodes, 1));
    for i = 2:length(p)
        if p(i) - p(i-1) < rtol
            id_toRemove(end+1) = topNodes(idx(i-1));
        end
    end

    % Remap T before removing nodes
    id_toRemove = unique(id_toRemove);
    keepNodes   = setdiff(1:size(X,1), id_toRemove);
    nodeMap     = zeros(size(X,1), 1);
    nodeMap(keepNodes) = 1:numel(keepNodes);

    % Remove elements that contain a removed node
    badElems = any(ismember(T, id_toRemove), 2);
    T = T(~badElems, :);
    T = nodeMap(T);   % renumber
    X = X(keepNodes, :);

    periodicMesh.X = X;
    periodicMesh.T = T;   % BAMG topology preserved, no delaunay
end