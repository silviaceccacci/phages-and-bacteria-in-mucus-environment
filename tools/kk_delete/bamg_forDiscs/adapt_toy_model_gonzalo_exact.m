function [Xh,Th,Uh] = adapt_toy_model_gonzalo_exact(X,T,u,data,model)

    global machine;
    switch lower(machine)
        case 'abel'
            bamg_exe = '!/usr/local/ff++/mpich3/bin/bamg';
        case 'gonzalo'
            %addpath('./bamg_interface/gonzalo')
            bamg_exe = '!bamg';
        otherwise
            error('adapt: not known machine')
    end
    
    mesh.X = X';
    mesh.T = T;
    
    switch model
        case 1
            fileName = './toy_model_gonzalo/exacta/exact1';
            fullname = [fileName '.msh'];
        case 2
            fileName = './toy_model_gonzalo/exacta/exact2B';
            fullname = [fileName '.msh'];
        case 3
            fileName = './toy_model_gonzalo/exacta/exact3';
            fullname = [fileName '.msh'];
        otherwise
            error('Case not implemented yet');
    end
    
    fileName = './temps/z_temp';
    
    msh    = printMSH(fileName,mesh);
    f_file = printSOL(fileName,u);
    
    %create the instruction to Bamg
    adaptativeStep = instructionBamgExact(bamg_exe,msh,data,f_file,fullname);
    
    % Per setejar tamany maxim i minim de malla:
    %   -hmin $HMIN -hmax 3
    % Quan calgui interpolar una solucio 
    %   usar -rbb i -wbb 
    % Altres coses:
    %   Podem posar geometria exacta:
    %    -g i el fitxer
    %    -errg 0.05
    %    La meva intencio inicial es no ser conforme amb el disc,
    %    xq es un model, no geometria. Aixi generalitza a 3D també
    
    fileName = './toy_model_gonzalo/exacta/exact_final';
    
    eval(adaptativeStep) %aquest funcio canvia directament els arxius z_temp.msh i z_temp.bb
    mesh = readBAMGmsh( fileName );
    
    Xh = mesh.X';
    Th = mesh.T;
    
    % read also the interpolated solution
    Uh = zeros(size(Xh,1),1); %aixo es canvia quan tinguem la solució
    
end

%     -b MESH_$j.msh -err 0.001 -errg 0.05 -AbsError \
% -hmin $HMIN -hmax 3 -Mbb MACH.bb -o MESH_$i.msh \
% -oamdba MESH_$i.amdba -ratio 2 -rbb SOL_$j.bb -wbb INIT_$i.bb

%adaptativeStep = ['!bamg -b ' msh ' -M ' mtr ' -o ' msh]
%adaptativeStep = ['!/usr/local/ff++/mpich3/bin/bamg -b ' msh ' -M ' mtr ' -o ' msh]

