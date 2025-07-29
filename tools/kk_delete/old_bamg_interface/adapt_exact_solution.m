% exact geometry
% model = 1 (linia {0}x[-25,25]), 2 (disc [-1,1]x[-25,25], 3 (2 discos [-40,-38]x[-25,25] i [40,42]x[-25,25])
% input: model
% output: Mesh (X,T)
function [Xh,Th,u,p,ux,uy] = adapt_exact_solution(...
    Xp,Tp,data,omega,adapt_variable,sol)

    error('BAMG does not allow adapt to sol with geometry.. need metric')

    [case_1disc,case_2disc,case_wods,case_wodsBig]=getAvailableCases();

    if omega.model==case_1disc
        
    else
        error('epp check adapt_exact')
    end
   
    %%
    [bamg_exe]=get_bamg_exe();

    
    %% geo mesh
    f_inp = './temps/exact_inp.msh';
    %f_out = './temps/exact_out.msh';
    
    hmin=data.hmin;
    hmax=data.hmax;
    
    x_LL = omega.x_LL ;
    x_LR = omega.x_LR ;
    y_LD = omega.y_LD ;
    y_LU = omega.y_LU ;
    X = [ x_LL y_LD 
          x_LR y_LD 
          x_LR y_LU
          x_LL y_LU ];
    XD = omega.position_discs;
        
    X = [ X ; XD];
    
    T = [ 1 2 ; 2 3; 3 4; 4 1];
    T = [T ; T+4];
    
    hgeo = [str2double(hmax) .* ones(4,1)
            str2double(hmin) .* ones(4,1)];

    meshGeo.X = X';
    meshGeo.T = T;
    
    printGeoMSH(f_inp,meshGeo)%,hgeo)
    
    %% adaptive mesh gen
    geo_string = [' -g ' f_inp ];
%     adaptativeStep = [ bamg_exe ' -g ' f_inp  ' -nbv ' data.nbv ' -o ' f_out ];
%     adaptativeStep = [ adaptativeStep ' > ./temps/out_exact.txt' ];
% 
%     eval(adaptativeStep)
% 
%     mesh = readBAMGmsh( f_out );
%  
%     Xh = mesh.X';
%     Th = mesh.T;
    
    %%
    fileName_mesh_input = './temps/temp_mesh_input';
    mesh.X = Xp';
    mesh.T = Tp;
    fn_mesh_inp  = printMSH(fileName_mesh_input,mesh);
    
    fn_mesh_out0 = './temps/temp_mesh_output';
    fn_mesh_out  = [fn_mesh_out0 '.msh'];

    fn_sol_inp      = './temps/temp_multiSol_inpput.bb';
    fn_sol_out      = './temps/temp_multiSol_output.bb';
    
    if adapt_variable.u==1 && adapt_variable.d ==1 && adapt_variable.p ==0
        f_coupez = Fun(X(:,1),X(:,2),data.beta,omega);
        target_sol = [sol.u f_coupez];
    elseif adapt_variable.u==1 && adapt_variable.d ==1 && adapt_variable.p ==1
        f_coupez = Fun(X(:,1),X(:,2),data.beta,omega);
        target_sol = [sol.u sol.p f_coupez];
    elseif adapt_variable.u==1 && adapt_variable.d ==0 && adapt_variable.p ==0
        target_sol = [sol.u];
    elseif adapt_variable.u==1 && adapt_variable.d ==0 && adapt_variable.p ==1
        target_sol = [sol.u sol.p];
    else
        error('not implemented strategy')
    end 
    printMultipleSOL(fn_sol_inp,target_sol);
    
    [adaptStep]=createBamgString(bamg_exe,fn_mesh_inp,fn_mesh_out,...
        fn_sol_inp,data,fn_sol_out,geo_string);

    %adaptStep = [adaptStep ' -g ' f_inp ];
    %adaptStep = [adaptStep ' > ./temps/out_bamg_adapt.txt'];%PETA NO SE XQ
    disp(adaptStep)
    eval(adaptStep)
    
    mesh = readBAMGmsh(fn_mesh_out0);
    
    Xh = mesh.X';
    Th = mesh.T;
    
    
    
    %%
   %error('check files')
    
    delete(fn_mesh_inp);
    delete(fn_mesh_out);
    %delete(fn_sol_out);
    
    u=[];
    p=[];
    ux=[]; 
    uy=[];
    
    delete(f_inp);
    %delete(f_out);
end