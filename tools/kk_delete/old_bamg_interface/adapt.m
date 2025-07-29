function [Xh,Th,Uh] = adapt(X,T,u)

    global machine;
    switch lower(machine)
        case 'abel'
            %bamg_exe = '!/usr/local/ff++/mpich3/bin/bamg';
            bamg_exe = '!/Applications/FreeFem++.app/Contents/ff-4.11/bin/bamg';
        otherwise
            error('adapt: not known machine')
    end 

    fileName = './temps/z_temp';
    
    mesh.X = X';
    mesh.T = T;
    
    msh    = printMSH(fileName,mesh);
    f_file = printSOL(fileName,u);

    adaptativeStep = [ ...
        bamg_exe ' -aniso ' ' -b ' msh ...
        ' -err 0.1 -AbsError ' ...
        '-Mbb ' f_file ...
        ' -o '  msh ]
%         ' -power 2 '...

       % '-hmin 1.0 -hmax 10.0'...
%     adaptativeStep = [ ...
%         bamg_exe ' -iso ' ' -b ' msh ...
%         ' -err 0.001 -AbsError ' ...
%         '-Mbb ' f_file ...
%         ' -o '  msh ]
    
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
    
    eval(adaptativeStep)
    
    mesh = readBAMGmsh( fileName );
    
    Xh = mesh.X';
    Th = mesh.T;
    
    % read also the interpolated solution
    Uh = zeros(size(Xh,1),1);
    
end

%     -b MESH_$j.msh -err 0.001 -errg 0.05 -AbsError \
% -hmin $HMIN -hmax 3 -Mbb MACH.bb -o MESH_$i.msh \
% -oamdba MESH_$i.amdba -ratio 2 -rbb SOL_$j.bb -wbb INIT_$i.bb

%adaptativeStep = ['!bamg -b ' msh ' -M ' mtr ' -o ' msh]
%adaptativeStep = ['!/usr/local/ff++/mpich3/bin/bamg -b ' msh ' -M ' mtr ' -o ' msh]

