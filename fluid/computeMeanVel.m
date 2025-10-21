function umean = computeMeanVel(mesh,domain,solution)

u = solution.u;

%u_scalar = sqrt(sum(u.*u,2));
u_scalar = u(:,1); % we want left to right vel

umean = mean(u_scalar);

