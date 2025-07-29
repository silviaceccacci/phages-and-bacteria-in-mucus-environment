function [Xh,Th] = adapt_nonexact_before_adaptation(X,T,u,data)
warning('---- AQUI TAMB HAURIA DE TRIAR ENTRE EL MEU AMB COMPLEXITAT I EL BAMG ----')

    [bamg_exe]=get_bamg_exe(); 

    fileName_mesh_input = './temps/temp_mesh_input';
    fileName_mesh_output = './temps/temp_mesh_output';
    fileName_metric_input = './temps/temp_metric_input';
    
    mesh.X = X';
    mesh.T = T;
    
    mesh_input = printMSH(fileName_mesh_input,mesh);
    mesh_output = './temps/temp_mesh_output.msh';
    metric_input = printSOL(fileName_metric_input,u);
    
    %create the instruction to Bamg
    adaptativeStep = instructionBamg(bamg_exe,mesh_input,mesh_output,metric_input,data);
    
    adaptativeStep=[adaptativeStep '> ./temps/out_mesh_creation_bamg.txt'];
    eval(adaptativeStep) %aquest funcio canvia directament els arxius z_temp.msh i z_temp.bb
    mesh = readBAMGmsh(fileName_mesh_output);
    
    Xh = mesh.X';
    Th = mesh.T;
    
    %remove the files
    delete(mesh_input);
    delete(mesh_output);
    delete(metric_input);
    
end

%     -b MESH_$j.msh -err 0.001 -errg 0.05 -AbsError \
% -hmin $HMIN -hmax 3 -Mbb MACH.bb -o MESH_$i.msh \
% -oamdba MESH_$i.amdba -ratio 2 -rbb SOL_$j.bb -wbb INIT_$i.bb

%adaptativeStep = ['!bamg -b ' msh ' -M ' mtr ' -o ' msh]
%adaptativeStep = ['!/usr/local/ff++/mpich3/bin/bamg -b ' msh ' -M ' mtr ' -o ' msh]

