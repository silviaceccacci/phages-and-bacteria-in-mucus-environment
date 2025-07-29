function [X,T,nodesDISC,elementsDISC,n,sol_interp] = adaptMeshtoSolution_complexity(...
    ifplot,exacta,data,X,T,Xp,Tp,adapt_variable,ifpreviousmesh,omega,mesh_solution)
%% debugging
out_debug_sizes = true;
%% warnings and errors
n = 0;  
if(exacta == 1) error('Not implemented'); end
if(ifpreviousmesh == 0) error('should not enter here...'); end
if(adapt_variable.comb == 1) error('Removed option, see adaptMesh_sec to get it'); end
        
%     ux_interp = [];
%     uy_interp = [];
%
%     if(exacta==1)
%         [X,T,u_inter,p_inter,u_interx,u_intery] =...
%         adapt_exact_solution(Xp,Tp,data,omega,adapt_variable,mesh_solution);
%     else
%         [X,T,u_inter,p_inter,u_interx,u_intery] =...
%         adapt_multipleSol(Xp,Tp,data,omega,adapt_variable,mesh_solution);
% 
%     end
%%
mesh_solution
mesh   = mesh_solution;

h_max = str2double(  data.hmax);
h_min = str2double(  data.hmin);
targetNodes = str2num(data.nbv);
%%
if adapt_variable.u==1 && adapt_variable.d ==1 && adapt_variable.p ==0
    f_coupez = Fun(Xp(:,1),Xp(:,2),data.beta,omega);
    u = mesh_solution.u;
    target_sol = [u f_coupez];
else
    error('not implemented in adaptMeshtoSolution_complexity')
end
h = computeMetric_fromScalars(mesh,target_sol);

if(out_debug_sizes)
    opt_parav.exportName = 'debugging';
    opt_parav.counter = 0;
    opt_parav.u = [mesh_solution.ux mesh_solution.uy];
    opt_parav.p = mesh_solution.p;
    opt_parav.f = Fun(Xp(:,1),Xp(:,2),omega.beta,omega);
    opt_parav.h = h;
    %opt_parav.d = elementsDISC;
    exportMeshParaview(Xp,Tp,opt_parav);
end
%% GRADATE METRIC
warning('compute metric gradation!!!')
%% Compute metric with the desired complexity
h = computeNewMetric_withDesiredComplexity(mesh,h,targetNodes);

C = computeMetric_complexity(mesh,h);
err_numNodTheoretical = norm(2*C-targetNodes);
if(err_numNodTheoretical>1e-12)
    disp(['The complexity of the new metric is:          ' num2str(C)])
    error('Unexpected complexity...')
end
disp(['The expected num nodes of the new metric are: ' num2str(2*C)])

if(out_debug_sizes)
    opt_parav.exportName = 'debugging';
    opt_parav.counter = opt_parav.counter+1;
    opt_parav.h = h;
    exportMeshParaview(Xp,Tp,opt_parav);
end
%% Bound metric max/min values
h = min(h,h_max);
h = max(h,h_min);
h = computeNewMetric_withDesiredComplexity(mesh,h,targetNodes);

% maybe I want to bound first the eigenvalues and then add complexity

if(out_debug_sizes)
    opt_parav.exportName = 'debugging';
    opt_parav.counter = opt_parav.counter+1;
    opt_parav.h = h;
    exportMeshParaview(Xp,Tp,opt_parav);
end
%% Things in mind
% 1. compute first and second order derivatives of the field
% 2. compute metric
% 3. compute size field (optional, but use it now to simplify)
% 4. compute complexity
% 5. reescale field to have desired complexity

% Things to think about and try:
%  - intersect target mesh sizes before or after complexity...
%    A: if I do it after I ensure both are reproduces
%    B: if I intsersect firts, then none is overresolved 
%  I would choose B at first, if it works fine I stay with it

%%

error(' call bamg with metric file')
% 1. build interface with bamg when size field is available
%    (it is easy since requires few options)

%%
sol_interp.u  = u_inter;
sol_interp.p  = p_inter;

if(isempty(ux_interp)==0)
    sol_interp.ux = u_interx;
    sol_interp.uy = u_intery;
    delete(ux_interp);
    delete(uy_interp);
end

nodesDISC    = FindNodesDisc(X,omega);
[elementsDISC,nodesDISC] = computeElementsDisc(X,T,omega,nodesDISC);

fprintf('  numElem: %d   numNodes: %d \n',size(T,1),size(X,1))


end