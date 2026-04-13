function periodicMesh = makePeriodicMeshStokes(mesh, hmin)
    % Enforce periodicity for a mesh on a rectangle with dimensions Lx x Ly.
    % Hole preservation: Delaunay + centroid-in-original-mesh filter.
    % Periodicity guarantee: boundary nodes are built symmetrically AND
    % are never dropped by the unused-node stripping step.

    X = mesh.X;
    T = mesh.T;

    x0 = min(X(:,1));  y0 = min(X(:,2));
    x1 = max(X(:,1));  y1 = max(X(:,2));
    Lx = x1 - x0;
    Ly = y1 - y0;

    X(:,1) = X(:,1) - x0;
    X(:,2) = X(:,2) - y0;

    Xorig   = X;
    TR_orig = triangulation(T, Xorig);

    tol  = 1e-10;
    rtol = hmin * 0.5;

    % ------------------------------------------------------------------ %
    % 1.  Build symmetric boundary node sets                               %
    % ------------------------------------------------------------------ %
    leftNodes   = find(abs(X(:,1))       < tol);
    rightNodes  = find(abs(X(:,1) - Lx)  < tol);
    bottomNodes = find(abs(X(:,2))       < tol);
    topNodes    = find(abs(X(:,2) - Ly)  < tol);

    yLR = sort(unique([X(leftNodes,2);  X(rightNodes,2)]));
    yLR = collapseClose(yLR, rtol);

    xBT = sort(unique([X(bottomNodes,1); X(topNodes,1)]));
    xBT = collapseClose(xBT, rtol);

    % ------------------------------------------------------------------ %
    % 2.  Replace all boundary nodes with the symmetric sets               %
    % ------------------------------------------------------------------ %
    allBoundary      = unique([leftNodes; rightNodes; bottomNodes; topNodes]);
    X(allBoundary,:) = [];

    newLeft   = [zeros(length(yLR),1),     yLR(:)               ];
    newRight  = [Lx*ones(length(yLR),1),   yLR(:)               ];
    newBottom = [xBT(:),                   zeros(length(xBT),1)  ];
    newTop    = [xBT(:),                   Ly*ones(length(xBT),1)];

    newBdry = unique([newLeft; newRight; newBottom; newTop], 'rows');
    X = [X; newBdry];

    % Record which indices are boundary nodes (they must never be dropped)
    nX           = size(X, 1);
    nInner       = nX - size(newBdry, 1);
    bdryIdx      = (nInner+1 : nX)';   % boundary nodes are appended at the end

    % ------------------------------------------------------------------ %
    % 3.  Delaunay + centroid filter                                       %
    % ------------------------------------------------------------------ %
    DT    = delaunayTriangulation(X);
    T_new = DT.ConnectivityList;

    cx = (X(T_new(:,1),1) + X(T_new(:,2),1) + X(T_new(:,3),1)) / 3;
    cy = (X(T_new(:,1),2) + X(T_new(:,2),2) + X(T_new(:,3),2)) / 3;

    tid    = pointLocation(TR_orig, cx, cy);
    inside = ~isnan(tid);
    T_new  = T_new(inside, :);

    % ------------------------------------------------------------------ %
    % 4.  Strip unused nodes BUT always keep boundary nodes                %
    % ------------------------------------------------------------------ %
    used           = false(nX, 1);
    used(T_new(:)) = true;
    used(bdryIdx)  = true;   % <-- this is the key fix

    map            = cumsum(used);
    periodicMesh.X = X(used, :);
    periodicMesh.T = map(T_new);
end


function v = collapseClose(v, rtol)
    changed = true;
    while changed
        changed = false;
        i = 1;
        while i < length(v)
            if v(i+1) - v(i) < rtol
                v(i)   = (v(i) + v(i+1)) / 2;
                v(i+1) = [];
                changed = true;
            else
                i = i + 1;
            end
        end
    end
end