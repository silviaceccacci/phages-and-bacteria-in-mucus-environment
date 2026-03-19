In this folder there are three main files which are needed to obtain the results of this work:

1) main_im2mesh.m 
- It takes the mucus image mucus1.png in the /input folder as an input and it generates the mesh used for FE fluid solver, mucus1_mesh.mat (output).
- Full im2mesh code in main_im2mesh.m
- If adaptation is desired (do_adapt_mesh set to true), the user must install freeFem to adapt the mesh to the image using BAMG, and add his installation path to get_bamg_exe.m. BAMG is not linked, but called externally if desired. Other mesh generation tools as Triangle, or BAMG can be used if mesh adaptation is desired. The adaptation step can be skipped (do_adapt_mesh=false) and the image itself is used then as a structured mesh. Adaptation is recommended to obtain a lower node count mesh.

2) main_flowSolver.m 
- It takes the generated mesh mucus1_mesh.mat and it computes the flow solution of the Stokes-Brinkman equation for a given Darcy number, to be specified in the variable parameters.minDarcyNum. For the Darcy number Da=1e-4, it produces the interpolant mucus1_uInterp_1e-4.mat (output). 
- It uses functions defined in the folder /fluid and /interpolation.

3) main_particleFromImage.m 
- It takes the interpolant mucus1_uInterp_1e-4.mat as an input and it computes the particles dynamics based on the Langevin framework, for a given number of phages and bacteria, to be specified as num_phages and num_bacteria. It produces the cumulative number of phages attached over time in the file n_phages_vs_time_vec.txt (output). 
- It uses functions defined in the folder /particles.

All the outputs are saved in the folder /output.

To run the code, first run main_im2mesh.m, then main_flowSolver.m and finally main_particleFromImage.m.