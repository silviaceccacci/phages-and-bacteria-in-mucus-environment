function [f_out]=interpolateMeshes(mesh_in,mesh_out,f_in)
%fprintf('  Interp solution...'); %tic;

fn_in = './interpolation/temps/temp_in.msh';
printMSH(fn_in,mesh_in);

fn_out = './interpolation/temps/temp_out.msh';
printMSH(fn_out,mesh_out);

[bamg_exe]=get_bamg_exe();

fn_sol_in = './interpolation/temps/temp_sol.bb';
fn_sol_out = './interpolation/temps/temp_sol_interp.bb';

printMultipleSOL(fn_sol_in,f_in);

adaptativeStep = [ bamg_exe ' -b ' fn_in ' -r ' fn_out  ...
    ' -err 0.00000000001 ' ...
    ' -rbb ' fn_sol_in ' -wbb ' fn_sol_out];

adaptativeStep = [ adaptativeStep ' > ./interpolation/temps/out_bamg_interp.txt'];

eval(adaptativeStep)

f_out = readBAMGbb(fn_sol_out)';

delete(fn_in)
delete(fn_out)
delete(fn_sol_in)
delete(fn_sol_out)