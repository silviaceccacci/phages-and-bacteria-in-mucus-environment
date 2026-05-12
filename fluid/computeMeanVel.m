function umean = computeMeanVel(mesh, domain, solution)

    X = mesh.X;
    T = mesh.T;
    u = solution.u;
    ux = u(:,1);  % x-velocity (left to right)

    % Area of each triangle
    x1 = X(T(:,1),1); y1 = X(T(:,1),2);
    x2 = X(T(:,2),1); y2 = X(T(:,2),2);
    x3 = X(T(:,3),1); y3 = X(T(:,3),2);

    areas = 0.5 * abs((x2-x1).*(y3-y1) - (x3-x1).*(y2-y1));

    % Average ux at each triangle centroid
    ux_elem = (ux(T(:,1)) + ux(T(:,2)) + ux(T(:,3))) / 3;

    % Total domain area from domain extents (already known)
    A_total = ((domain.x_LR - domain.x_LL) * (domain.y_LU - domain.y_LD))*10^(-12);

    % Area-weighted integral divided by total area
    umean = sum(areas .* ux_elem) / A_total;

end