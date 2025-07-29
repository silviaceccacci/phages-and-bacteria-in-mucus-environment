function [mesh3D] = generateTetMesh(surfaceMesh,parameters,tempFolder)

    tempName = [tempFolder 'ztet_temp'];
        
    [inpName] = printMeshToTetGen(surfaceMesh,tempName,parameters);

    logName = [tempName '.log'];
    
    tetExe = '!./src/tetMeshing/tetgen/tetgen ';
     
    % ==> frmo Discmesh with surface constrain
    %stringSize = '100';
    %optionsTetgen = [' -p -O9/7 -q1.2/15 -a' stringSize ' -Y -V -BF '  inpName '>' logName ];
    
    % ==> basic Delaunay
    %optionsTetgen = [' -p -V -BF '  inpName '>' logName ];
    
    faceSwitch = ' -f ';
    
    optionsTetgen = [' -p -q1.2 -V -BF ' inpName ];%'>' logName ];
    optionsTetgen = [' -p -V -BF ' inpName ];%'>' logName ];
    optionsTetgen = [' -p -q2/10 -V -BF ' inpName ];%'>' logName ];
    optionsTetgen = [' -p -q1.2/15 -V -BF ' inpName ];%'>' logName ];
    optionsTetgen = [' -p  -O9/7 -q2/10 -V -m -nn ' inpName ];%'>' logName ];
    optionsTetgen = [' -p  -O9/7 -q1.2/15 -V -m -nn ' inpName ];%'>' logName ];
    optionsTetgen = [' -p  -O9/7 -q2/10 -V -m -nn ' inpName ];%'>' logName ];
    optionsTetgen = [' -p  -O9/7 -q2/10 -V -m -nn ' faceSwitch inpName ];%'>' logName ];
    optionsTetgen = [' -p  -O9/7 -q2/10 -m -nn ' faceSwitch inpName ];%'>' logName ];
%     optionsTetgen = [' -p -d -O9/7 -q2/10 -m -nn ' faceSwitch inpName ];%'>' logName ];
    %optionsTetgen = [' -p  -V -nn ' inpName ];%'>' logName ];
    %optionsTetgen = [' -p  -O9/7 -q -V -m -nn ' inpName ];%'>' logName ];
    %optionsTetgen = [' -pqm -V -nn ' inpName ];%'>' logName ];
    
%     stringSize = '100';
%     optionsTetgen = [' -p -O9/7 -q1.2/15 -a' stringSize ' -V -BF '  inpName '>' logName ];    
    
%     warning('HEI, I am changing the mesher to do CDT without quality')
%     optionsTetgen = [' -p  -m -nn ' faceSwitch inpName ];%'>' logName ];

    tetMeshing = [ tetExe optionsTetgen ]% '>' logName ] 
    
    eval(tetMeshing)
    
    [mesh3D]=readMeshFromTetGen(tempName,faceSwitch);
    
    mesh3D.name = [surfaceMesh.name '_tet3D'];
end


% AIXO ES DE DISCMESH.. PERO ARA NO HAIG DE FICAR CONSTRAIN AL CONTORN,
% SINO QUE LI PUC DEIXAR SER FELISSSS
%     call system(TRIM(tetMesherExe)//&
%       ' -p -O9/7 -q1.2/15 -a'//TRIM(stringTetSize)//' -Y -V -BF '//gapFile&
%       // ' &> '//TRIM(path)//TRIM(name)//'.'//TRIM(gapName)//'.Meshing.log' ) !! USE BETTER THIS ONE

