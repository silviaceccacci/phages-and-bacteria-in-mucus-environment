function [Xh,Th,u,p,ux,uy] = adapt_multipleSol(X,T,data,omega,adapt_variable,sol)

    [bamg_exe]=get_bamg_exe();

    fileName_mesh_input = './temps/temp_mesh_input';
    mesh.X = X';
    mesh.T = T;
    fn_mesh_inp  = printMSH(fileName_mesh_input,mesh);
    
    fn_mesh_out0 = './temps/temp_mesh_output';
    fn_mesh_out  = [fn_mesh_out0 '.msh'];

    fn_sol_inp      = './temps/temp_multiSol_inpput.bb';
    fn_sol_out      = './temps/temp_multiSol_output.bb';
    
    if adapt_variable.u==1 && adapt_variable.d ==1 && adapt_variable.p ==0
        f_coupez = Fun(X(:,1),X(:,2),data.beta,omega);
        u = sol.u;
        %u = u.^4;
        %u = tanh(u*10)
        target_sol = [u f_coupez];
        %target_sol = sol.u.*sqrt(f_coupez);
    elseif adapt_variable.u==1 && adapt_variable.d ==1 && adapt_variable.p ==1
        f_coupez = Fun(X(:,1),X(:,2),data.beta,omega);
        target_sol = [sol.u sol.p f_coupez];
    elseif adapt_variable.u==1 && adapt_variable.d ==0 && adapt_variable.p ==0
        target_sol = [sol.u];
        %target_sol = [sol.ux sol.uy];
%         f_u = tanh(abs(sol.u-0.8));
%         f_coupez = Fun(X(:,1),X(:,2),data.beta,omega);
%         a = 0.9;
%         b = 0.1;
%         %f_u = tanh(1.5*abs(sol.u-0.0));
%         %a = 1;
%         %b = 0;
%         target_sol = [a.*f_u + b.*f_coupez];
%         disp('[a.*f_u + b.*f_coupez]')
    elseif adapt_variable.u==1 && adapt_variable.d ==0 && adapt_variable.p ==1
        target_sol = [sol.u sol.p];
    else
        error('not implemented strategy')
    end 
    printMultipleSOL(fn_sol_inp,target_sol);
    
    [adaptativeStep]=createBamgString(bamg_exe,fn_mesh_inp,fn_mesh_out,...
        fn_sol_inp,data,fn_sol_out);

    adaptativeStep = [adaptativeStep ' > ./temps/out_bamg_adapt.txt'];%PETA NO SE XQ
    %disp(adaptativeStep)
    eval(adaptativeStep)
    
    mesh = readBAMGmsh(fn_mesh_out0);
    
    Xh = mesh.X';
    Th = mesh.T;
    
    delete(fn_mesh_inp);
    delete(fn_mesh_out);
    %delete(fn_sol_out); % not created right now
    
    u=[];
    p=[];
    ux=[]; 
    uy=[];
    
    % HERE I CAN READ INTERP SOLUTIONS
%     u  = readBAMGbb(fn_sol_output);
%     p  = readBAMGbb(fn_press_output);
%     
%     delete(mesh_input);
%     delete(mesh_output);
%     delete(solution_output);
%     delete(pressure_output);
    
end