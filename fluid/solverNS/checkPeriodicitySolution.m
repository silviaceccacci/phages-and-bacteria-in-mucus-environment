function isPeriodic=checkPeriodicitySolution(sol,list_nodes,map_to_master)
isPeriodic = true;
for inode = 1:length(list_nodes)

    slave = list_nodes(inode);
    master = map_to_master(slave);
    
    sol_1 = sol(slave,:);
    sol_2 = sol(master,:);

    e_per = norm(sol_1-sol_2);
    if(e_per>1e-10)
%         sol_1
%         sol_2
%         e_per
        isPeriodic =false;
    end
end

