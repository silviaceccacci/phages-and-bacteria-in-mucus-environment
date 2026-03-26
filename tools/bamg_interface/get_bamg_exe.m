function [bamg_exe]=get_bamg_exe()

global machine;
switch lower(machine)
    case 'abelm2'
        bamg_exe = '!/Applications/FreeFem++.app/Contents/ff-4.11/bin/bamg';
    case 'abel'
        bamg_exe = '!/usr/local/ff++/mpich3/bin/bamg';
    case 'silvia'
        bamg_exe = '!/usr/local/ff++/mpich3/bin/bamg';
    otherwise
        error('get_bamg_exe: not known machine')
end
