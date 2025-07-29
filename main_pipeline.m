clear all; close all;

disp('Zero step: ????')
disp(' chose image')
disp(' filter image')

disp('First step: ')
disp('   Generate mesh from image')
disp('   Set mesh parameters in main_im2mesh')
main_im2mesh

disp('Second step:')
disp('   Run flow solver to compute a flow velocity field')
main_flowSolver

disp('Third step:')
disp('   Compute particle dynamics')
main_particleFromImage

disp('____________________________________')
disp('Possible changes: ')
disp('- Particles could influence flow, so we could loop flow-particles')
disp('- We could then also adapt the mesh to the particles too to capture them')