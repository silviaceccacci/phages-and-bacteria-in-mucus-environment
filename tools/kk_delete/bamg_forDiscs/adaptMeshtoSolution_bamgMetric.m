function [X,T,nodesDISC,elementsDISC,n,sol_interp] = adaptMeshtoSolution_bamgMetric(...
    ifplot,exacta,data,X,T,Xp,Tp,adapt_variable,ifpreviousmesh,omega,mesh_solution)
    n = 0;  
    %if(exacta == 1) error('Not implemented'); end
    if(ifpreviousmesh == 0) error('should not enter here...'); end
    if(adapt_variable.comb == 1) error('Removed option, see adaptMesh_sec to get it'); end
        
%     fileName_metric_input      = './temps/temp_metric_input';
%     fileName_solution_input    = './temps/temp_solution_input';
%     fileName_pressure_input    = './temps/temp_pressure_input';
%     fn_ux_inp    = './temps/temp_ux_inp';
%     fn_uy_inp    = './temps/temp_uy_inp';
%    
%     combination_input = ' ';
%     solution_input = printSOL(fileName_solution_input,mesh_solution.u');
%     pressure_input = printSOL(fileName_pressure_input,mesh_solution.p');
%     U2 = Fun(X(:,1),X(:,2),data.beta,omega);
%     metric_input = printSOL(fileName_metric_input,U2');
%     ux_interp = [];%printSOL(fn_ux_inp,mesh_solution.ux');
%     uy_interp = [];%printSOL(fn_uy_inp,mesh_solution.uy');
% 
%     [X,T,u_inter,p_inter,u_interx,u_intery] =...
%         adapt_nonexact_after_adaptation(Xp,Tp,data,...
%         adapt_variable,metric_input,solution_input,...
%         pressure_input,combination_input,...
%         ux_interp,uy_interp);
    
    ux_interp = [];
    uy_interp = [];

    if(exacta==1)
        error('not available option..')
        [X,T,u_inter,p_inter,u_interx,u_intery] =...
        adapt_exact_solution(Xp,Tp,data,omega,adapt_variable,mesh_solution);
    else
        [X,T,u_inter,p_inter,u_interx,u_intery] =...
        adapt_multipleSol(Xp,Tp,data,omega,adapt_variable,mesh_solution);

    end
    
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