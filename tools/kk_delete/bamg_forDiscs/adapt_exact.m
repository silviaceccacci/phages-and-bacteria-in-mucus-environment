% exact geometry
% model = 1 (linia {0}x[-25,25]), 2 (disc [-1,1]x[-25,25], 3 (2 discos [-40,-38]x[-25,25] i [40,42]x[-25,25])
% input: model
% output: Mesh (X,T)
function [Xh,Th] = adapt_exact(omega,data)

    [case_1disc,case_2disc,case_wods,case_wodsBig]=getAvailableCases();

    if omega.model==case_1disc
        
    else
        error('epp check adapt_exact')
    end
   
    %%

    [bamg_exe]=get_bamg_exe();

    f_inp = './temps/exact_inp.msh';
    f_out = './temps/exact_out.msh';
    
    %% geo mesh
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
    
%     XD2= XD;
%     l = 2;
%     XD2(1,:) = XD2(1,:) - [l,l];
%     XD2(2,:) = XD2(1,:) + [l,-l];
%     XD2(3,:) = XD2(1,:) + [l,l];
%     XD2(4,:) = XD2(1,:) + [-l,l];
%    
%     X = [ X ; XD; XD2];
%     
%     T = [ 1 2 ; 2 3; 3 4; 4 1];
%     T = [T ; T+4; T+8];
%     
%     hgeo = [str2num(hmax) .* ones(4,1)
%             str2num(hmin) .* ones(4,1)
%             str2num(hmax) .* ones(4,1)];
    
    %X = [X; XD2];
    %hgeo = [hgeo; str2num(hmax) .* ones(4,1)];

    meshGeo.X = X';
    meshGeo.T = T;
    
    printGeoMSH(f_inp,meshGeo,hgeo)
    
    %% adaptive mesh gen
   
    adaptativeStep = [ bamg_exe ' -g ' f_inp  ' -nbv ' data.nbv ' -o ' f_out ];
    adaptativeStep = [ adaptativeStep ' > ./temps/out_exact.txt' ];

    eval(adaptativeStep)

    mesh = readBAMGmsh( f_out );
 
    Xh = mesh.X';
    Th = mesh.T;
    
    
    %%
%    error('check files')
    
    delete(f_inp);
    delete(f_out);
end