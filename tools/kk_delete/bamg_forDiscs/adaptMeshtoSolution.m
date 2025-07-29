function [X,T,nodesDISC,elementsDISC,n,sol_interp] = adaptMeshtoSolution(...
ifplot,exacta,data,X,T,Xp,Tp,adapt_variable,ifpreviousmesh,omega,mesh_solution)

metric_bamg = 1;
metric_complexity = 2;

metricConstructionType = metric_bamg;

switch metricConstructionType
    case metric_bamg
        [X,T,nodesDISC,elementsDISC,n,sol_interp] = adaptMeshtoSolution_bamgMetric(...
        ifplot,exacta,data,X,T,Xp,Tp,adapt_variable,ifpreviousmesh,omega,mesh_solution);
    case metric_complexity
        [X,T,nodesDISC,elementsDISC,n,sol_interp] = adaptMeshtoSolution_complexity(...
        ifplot,exacta,data,X,T,Xp,Tp,adapt_variable,ifpreviousmesh,omega,mesh_solution);
    otherwise
        error('not implemented strategy')
end


end