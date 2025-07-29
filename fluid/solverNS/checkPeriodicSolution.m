function []=checkPeriodicSolution(mesh,solution,BC)

isPer_uTB = true; isPer_uLR=true; isPer_pTB=true; isPer_pLR = true;
if(BC.periodic.topBot)
    isPer_uTB=checkPeriodicitySolution(solution.u ,mesh.periodic_maps.top ,mesh.periodic_maps.topBottom);
end
if(BC.periodic.leftRight)
    isPer_uLR=checkPeriodicitySolution(solution.u ,mesh.periodic_maps.left,mesh.periodic_maps.leftRight);
end
if(BC.periodic.p.topBot)
    isPer_pTB=checkPeriodicitySolution(solution.p0,mesh.periodic_maps.top ,mesh.periodic_maps.topBottom);
end
if(BC.periodic.p.leftRight)isPer_pLR=checkPeriodicitySolution(solution.p0,mesh.periodic_maps.left,mesh.periodic_maps.leftRight);
    isPer_pLR=checkPeriodicitySolution(solution.p,mesh.periodic_maps.top ,mesh.periodic_maps.topBottom);
end
isPerSol = [isPer_uTB, isPer_pTB, isPer_uLR, isPer_pLR];
if(any(isPerSol==false))
    if(isPer_uTB==false || isPer_uLR==false), disp('Not periodic vel'),end
    if(isPer_pTB==false || isPer_pLR==false), disp('Not periodic press'),end
    error('not periodiccccccc when it should be!!!!')
end
