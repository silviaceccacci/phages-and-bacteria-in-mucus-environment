% Solve Stokes / Navier Stokes equations with permeability (NL with Picard)
% MINI for the mesh discretization vel/pers
%
function [solution] = NavierStokesMini(mesh,parameters)
if(nargin<2)
    error('call NavierStokesMini with 2 inputs!!!')
end
%% PARAMETERS
TOL              = parameters.TOL_solver;
model            = parameters.model;
do_convective    = parameters.convective;
do_parallAssembly= parameters.is_parall;
do_dirichletLeft     = parameters.BC.dirichletLeft;
do_dirichletTopBot   = parameters.BC.dirichletTopBot;
do_slidingTopBot     = parameters.BC.slidingTopBot;
isPeriodic_topBot    = parameters.BC.periodic.topBot;
isPeriodic_leftRight = parameters.BC.periodic.leftRight;
isPeriodic = isPeriodic_leftRight || isPeriodic_topBot;

boundaryValue    =@model.u_dir; 

if(do_convective)
    ifPicard = true;
    maxIterNonLinear = 50 ; %maximum iterations solver
else
    maxIterNonLinear = 1;
end
%% Mesh variables
domain = mesh.domain;
X =mesh.X;
T =mesh.T;
Xp=mesh.Xp;
Tp=mesh.Tp;
nOfNodes = size(X,1)+size(T,1) ;
nOfNodesp = size(Xp,1);

node_to_dof_vx = @(nodes)(nodes           );
node_to_dof_vy = @(nodes)(nodes+  nOfNodes);
node_to_dof_p  = @(nodes)(nodes+2*nOfNodes);
node_to_dofs_v = @(nodes)([node_to_dof_vx(nodes);node_to_dof_vy(nodes)]);
node_to_dofs   = @(nodes)([node_to_dofs_v(nodes);node_to_dof_p(nodes)]);
%% INITIAL CONDITION: u0
u0=zeros(nOfNodes,2); %u0(1:nOfNodesp,2) = readBAMGbb(fileName_ic_uy);

% Interpolate on the midle of the element for the MINI additional dof
u0((nOfNodesp+1):end,:) = ( u0(T(:,1),:)+u0(T(:,2),:)+u0(T(:,3),:) )/3.0;
%% Definition of boundary conditions 
x = X(:,1); y = X(:,2); 
tol=1.e-5; 
L_char = sqrt( (domain.x_LR-domain.x_LL)^2 + (domain.y_LU-domain.y_LD)^2 );
tol=tol*L_char;
if(do_dirichletLeft && do_dirichletTopBot)
    if(isPeriodic)
        error('Error considering all dirichlet with periodicity')
    end
    %Nodes on the boundary for Dirichlet BC (up,down,left) we have not put the right border (Neumann)
    nodesCCD = find(abs(x-domain.x_LL)<tol|abs(y-domain.y_LU)<tol|abs(y-domain.y_LD)<tol); 

    XnodesCCD  = X(nodesCCD,:);  %coordinates of the nodes on the boundary
    dofsCCD    = [nodesCCD; nodesCCD+nOfNodes];
    uCCD       = boundaryValue(XnodesCCD);%,domain);
    
    numNodesCCD = length(nodesCCD);
    u0(nodesCCD,1) = uCCD(1:numNodesCCD) ;
    u0(nodesCCD,2) = uCCD(numNodesCCD+(1:numNodesCCD)) ;

elseif(do_dirichletLeft && do_dirichletTopBot==false)
    if(do_slidingTopBot)
        % Dirichlet only on inflow, rest: symm 
        nodes_sliding = find(abs(x-domain.x_LL)>tol & ...
            (abs(y-domain.y_LU)<tol|abs(y-domain.y_LD)<tol));
    elseif(isPeriodic_topBot)
        nodes_sliding = zeros(0,1);
    else
        error('not possible')
    end

    uCCD_slide  = zeros(length(nodes_sliding),1);
    u0(nodes_sliding,2) =  uCCD_slide;

    nodes_inflow  = find(abs(x-domain.x_LL)<tol);
    dofsCCD     = [nodes_inflow; (nodes_inflow+nOfNodes); (nodes_sliding+nOfNodes)];
    uCCD_inflow = boundaryValue(X(nodes_inflow,:));
    numNod_inflow = length(nodes_inflow);
    u0(nodes_inflow,1)  = uCCD_inflow(1:numNod_inflow) ;
    u0(nodes_inflow,2)  = uCCD_inflow(numNod_inflow+(1:numNod_inflow)) ;
    uCCD = [uCCD_inflow ; uCCD_slide];

    nodesCCD = nodes_inflow;
    prescribedValues = uCCD;
    prescribedDOF    = dofsCCD;
elseif(do_dirichletLeft==false && do_dirichletTopBot==false)
    if isfield(parameters.BC,'noSlipWalls') && parameters.BC.noSlipWalls ...
            && isfield(mesh,'solid_wall_nodes')
        nodesCCD         = mesh.solid_wall_nodes;
        prescribedDOF    = [nodesCCD; nodesCCD + nOfNodes];
        prescribedValues = zeros(2*length(nodesCCD), 1);
        u0(nodesCCD, :)  = 0;
        disp('  No-slip BC applied on mucin walls')
    else
        nodesCCD         = [];
        prescribedValues = zeros(0,1);
        prescribedDOF    = zeros(0,1);
        disp('  No BC on mucin walls (free slip)')
    end
end
clear x; clear y;

dof_p_to_fix = [];
value_p_to_fix =[];
if(isPeriodic_leftRight && isPeriodic_topBot) 
    disp('Fixing a pressure dof')
    xcm = sum(Xp,1)/size(Xp,1);
    d   = [Xp(:,1)-xcm(1),Xp(:,2)-xcm(2)];
    d   = sum(d.*d,2);
    [mind,central_p_node]=min(d);
    nodeToFixP   = central_p_node;
    dof_p_to_fix = node_to_dof_p(nodeToFixP);
    value_p_to_fix = 0;
    clear xcm; clear d;
%     if(parameters.BC.periodic.p.lagrange)
%         error('check if i need to fix pressure dof')
%     end
end

periodic_data=compute_dofs_periodic(mesh,parameters);

model.periodic_data = periodic_data;
slaves_periodic=periodic_data.slaves_periodic;
master_periodic=periodic_data.master_periodic;
slaves_pjump=periodic_data.slaves_pjump;
master_pjump=periodic_data.master_pjump;
if(isempty(intersect(slaves_periodic,master_periodic))==false)
    intersect(slaves_periodic,master_periodic)
    error('slaves are masters!?!?!?!')
end
if(isempty(intersect(slaves_pjump,master_pjump))==false)
    intersect(slaves_periodic,master_periodic)
    error('slaves are master_pjump!?!?!?!')
end
%% System computation
fprintf('  CreMat... '); tic;
[K,G,f] = CreMat_par_faster(X,T,Xp,Tp,1,model,do_parallAssembly);
toc

isConvergedNonLinearTerm = false;
while(isConvergedNonLinearTerm==false)
    %% Picard -> update convective term
    if(do_convective)
        error('ensure convective terms is well implemented (periodic,etc)')
        disp('  CreConv...');  tic;
        C=CreConv_par_faster(u0,X,T,1,do_parallAssembly);
        toc;
    else
        C = sparse(size(K,1),size(K,2));
    end    
    %% Sytem matrix
    numConstrains = length(slaves_pjump);
    if(parameters.BC.periodic.p.lagrange)
        fprintf('  Using Lagrange multipliers to impose Delta p\n')
        fprintf('  Assemble [K+M+C G 0; G 0 L; 0 L 0]... '); tic;
        L = sparse(numConstrains,nOfNodesp); %size(G,1)  
        pjump = parameters.pressureJump*ones(numConstrains,1);  
        for islave = 1:length(slaves_pjump)
            theSlave  = slaves_pjump(islave)-2*nOfNodes; % dofs in pres
            theMaster = master_pjump(islave)-2*nOfNodes; % dofs in pres
            L(islave,theSlave) =  1.0;
            L(islave,theMaster) = -1.0;
        end
        myZerosGG = sparse(size(G,1),size(G,1));
        myZerosLG = sparse(size(L,1),size(G,2));
        myZerosLL = sparse(size(L,1),size(L,1));
        A = [   K+C         G'     myZerosLG'
                 G      myZerosGG     L'
             myZerosLG      L      myZerosLL];
        b = [        f
            zeros(nOfNodesp,1)
                   pjump       ];
        clear C; clear pjump; clear L; clear myZerosLL; clear myZerosLG;
    else
         fprintf('  Assemble [K+M+C G;G 0]... '); tic;
         myZeros = sparse(nOfNodesp,nOfNodesp);
         A= [K+C G';G myZeros];
         b = [f; zeros(nOfNodesp,1)];
         clear myZeros; clear C;
    end

    %% Solve
    toc; fprintf('  Set unkowns and dirichlet on source term... '); tic;
    %__Imposing of Dirichlet boundary conditions (system reduction)
    % and computing actual degrees of freedom (not boundary nodes, or periodic)   
    unknowns = 1:size(A,1); % 2*nOfNodes+nOfNodesp
    unknowns = setdiff(unknowns,prescribedDOF); 
    unknowns = setdiff(unknowns,slaves_periodic); %periodic_slave_dofs
    unknowns = setdiff(unknowns,dof_p_to_fix);
    prescribedDOF = [prescribedDOF; dof_p_to_fix];
    prescribedValues = [prescribedValues; value_p_to_fix];

    b = b(unknowns)-A(unknowns,prescribedDOF)*prescribedValues;
    A = A(unknowns,unknowns);
    toc;

%     tic; 
%     condA = condest(A);
%     fprintf('  Condition number: %e ...',condA); toc;
%     if(condA>1e15) 
%         error('Too large condition number.. maybe BCs are not correct? lack of p fixed? :-)')
%     end

    %__System solution
    fprintf('  Solve linear system... '); tic;
    perm      = symrcm(A); %Symmetric reverse Cuthill-McKee permutation
    sol       = zeros(size(A,1),1); 
    sol(perm) = A(perm,perm)\b(perm);
    clear A; toc;

    auxsol                = zeros(2*nOfNodes+nOfNodesp+numConstrains,1);
    auxsol(unknowns)      = sol                          ;
    auxsol(prescribedDOF) = prescribedValues             ;
    if(isPeriodic)
        disp('  Assigning periodic dof values...')
        auxsol(slaves_periodic) = auxsol(master_periodic);
    end
    auxsol(dof_p_to_fix)=value_p_to_fix;

    u = auxsol(1:2*nOfNodes);  
    u = reshape(u,nOfNodes,2);
    p = auxsol(2*nOfNodes+1:(2*nOfNodes+nOfNodesp));
    if(parameters.BC.periodic.p.lagrange)
        lambda = auxsol((2*nOfNodes+nOfNodesp+1):end);
    end
    
    if(do_convective)
        err1_u = max( sqrt(sum((u0-u).^2,2)) ) / max( sqrt(sum(u0.^2,2)) );
        err1_p = max(abs(p0-p))/max(abs(p(:)));
        maxE = max(err1_u,err1_p);
        
        u0 = u; p0 = p;
    
        fprintf('  Iteration %d -> Er: %f\n',i,maxE);
        
        isLowError = maxE<TOL;
        isIterViolated = iter>maxIterNonLinear;
        isConvergedNonLinearTerm = isLowError | isIterViolated;
    else
        isConvergedNonLinearTerm = true;
    end
end
clear K; clear G;

%if(isPeriodic_leftRight)
    imposedGrad = parameters.pressureGrad;
    imposedGrad = imposedGrad*Xp(:,1);
    p_solved = p;
    p = p+imposedGrad;
%end

%POSTPROCESS
u = u(1:nOfNodesp,:);
ux = u(:,1);
uy = u(:,2);

solution.u    = u;
solution.p    = p;
solution.umag = sqrt(ux.^2+uy.^2);
solution.p0 = p_solved;

node_marks = zeros(size(Xp,1),1);
if(do_slidingTopBot)
    node_marks(nodes_inflow) = 1;
    node_marks(nodes_sliding) = 3;
else
    node_marks(nodesCCD) = 1;
end
solution.node_marks = node_marks;

if(do_convective)
    if isLowError
        %fprintf('Convergence in %3d steps\n',i);
    else
        global case_name;
        opt_parav_debug.exportName = [case_name '_notConverged'];
        opt_parav_debug.u    = solution.umag;
        opt_parav_debug.p    = solution.p;
        exportMeshParaview(Xp,Tp,opt_parav_debug);
        error('No convergence in NavierStokesDiscMini!!!!!!!');
    end
end


end

%             error('not sym matrix')
%             fprintf('  Assemble [K+M+C G;G 0]... '); tic;
%             myZeros = sparse(nOfNodesp,nOfNodesp);
%             A= [K+C G';G myZeros];
%             b = [f; zeros(nOfNodesp,1)];
%             clear myZeros; clear C;
% 
%             for islave = 1:length(slaves_pjump)
%                 theSlave = slaves_pjump(islave);
%                 theMaster = master_pjump(islave);
%                 A(theSlave,:) = 0.0;
%                 A(theSlave,theSlave) =  1.0;
%                 A(theSlave,theMaster) = -1.0;
%                 b(theSlave) = parameters.pressureJump;
%             end

