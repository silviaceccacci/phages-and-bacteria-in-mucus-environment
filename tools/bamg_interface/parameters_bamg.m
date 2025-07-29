function [data] = parameters_bamg()
%% not interesting parameters
data.coef     = '1';      %default: 1
data.ratio    = '0';      %default: 0
data.errg     = '0.1';    %Geometry error, never bigger than 0.7071

%% error parameters
data.err      = '0.1';    %p1 interpolation error, default 0.01
data.iso      = '0';      %1 si volem que sigui isotropic, 0 vol dir anisotropic
data.anisomax = '100';     %0 without anisomax

%% mesh size parameters
data.nbv    = '1000000';

data.hmin   = '1.0'   ;
data.hmax   = '100.0'   ;
% data.hmin   = '0'   ;
% data.hmax   = '0'   ;


% -AbsError 
% -NoRescaling 
% -NbJacobi 2 
% -NbSmooth 5 
% -hmax 2 
% -hmin 0.0000005 
% -ratio 0 
% -nbv 100000 
% -v 4 
% -err 0.05 
% -errg 0.01