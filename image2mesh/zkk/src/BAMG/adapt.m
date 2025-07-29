function [Xh,Th] = adapt(X,T,h,options)

    if(nargin>3 && isfield(options,'fileName'))
        fileName = options.fileName;
    else
        fileName = 'zbamg_temp';
    end
    if(nargin>3 && isfield(options,'meshToAdapt'))
        meshToAdapt = options.meshToAdapt;
        twoMeshes = true;
    else
        twoMeshes = false;
    end

    mesh.X = X';
    mesh.T = T;
    
    msh = printMSH(fileName,mesh);
    mtr = printMTR(fileName,h);
    
    
    if(nargin>3 && isfield(options,'scalar') )
        scalarFileName = printFIELD(fileName,options.scalar);
    end
    
    %BAMGoptions = ' -ratio 50 -anisomax 1000000 -hmax 10 -hmin 1';
    
%     BAMGoptions = [' -AbsError -NoRescaling -NbJacobi 2 -NbSmooth 5 '...
%                    '-hmax 10 -hmin 0.1 -ratio 0 -nbv 100000 ' ...
%                    '-v 4 -err 0.05 -errg 0.01 '];  
               
% 	BAMGoptions = [' -AbsError -NoRescaling -NbJacobi 2 -NbSmooth 5 '...
%                    '-hmax 10 -hmin 0.1 -ratio 0 -nbv 100000 ' ...
%                    '-v 4 -err 0.05 -errg 0.01  -Mbb ' scalarFileName];
%                
% 	BAMGoptions = [' -RelError -NoRescaling -NbJacobi 5 -NbSmooth 10 '...
%                    ' -hmax 20 -ratio 50 -err 0.5 ' ...
%                    '-v 4  -Mbb ' scalarFileName];
               
% 	BAMGoptions = [' -AbsError -NoRescaling -NbJacobi 2 -NbSmooth 5 '...
%                    ' -hmax 10 -ratio 50 -v 4 -err 0.01 -Mbb ' scalarFileName];
	BAMGoptions = [' -iso -AbsError -NoRescaling -NbJacobi 2 -NbSmooth 5 '...
                   ' -power 4 -ratio 100 -nbv 1000000 ' ...
                   '-v 4 -err 0.05 -Mbb ' scalarFileName];
               
%                -AbsError -NoRescaling -NbJacobi 2 -NbSmooth 5 -hmax 2 -hmin 0.0000005 -ratio 0 -nbv 100000 -v 4 -err 0.05 -errg 0.01
    
    if(twoMeshes)
        error('Wrong option')
    else
        global machine;
        if(strcmp(machine,'abel'))
            BAMGroute = '!/usr/local/ff++/mpich/3.41/bin/bamg ';
        elseif(strcmp(machine,'xevi'))
            BAMGroute = '!/usr/local/ff++/mpich/3.43-1/bin/bamg ';
        else
            error('Not identified machine');
        end
    end
    
    %adaptativeStep = [ BAMGroute ' -b ' msh ' -M ' mtr  BAMGoptions  ' -o ' msh ] ;
    adaptativeStep = [ BAMGroute ' -b ' msh  BAMGoptions  ' -o ' msh ] ;
    
    eval(adaptativeStep)
    mesh = readBAMGmsh( fileName );
    
    Xh = mesh.X';
    Th = mesh.T;
    
end