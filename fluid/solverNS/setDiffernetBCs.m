function [BC_dirichletLeft,BC_lagrangeJump,BC_pgrad,BC_pgrad_noperx]=setDiffernetBCs()
% This works fine and gives a well conditioned matrix:
BC.dirichletLeft        = true;
BC.dirichletTopBot      = false;
BC.slidingTopBot        = false;
BC.periodic.topBot      = true;
BC.periodic.leftRight   = false;
BC.periodic.p.topBot    = true;
BC.periodic.p.leftRight = false;
BC.periodic.p.lagrange  = false;
BC.label = 'dir_per';
BC_dirichletLeft = BC;

% BC.dirichletLeft        = true;
% BC.dirichletTopBot      = true;
% BC.slidingTopBot        = false;
% BC.periodic.topBot      = false;
% BC.periodic.leftRight   = false;
% BC.periodic.p.topBot    = false;
% BC.periodic.p.leftRight = false;
% BC.periodic.p.lagrange  = false;
% BC.label = 'dir_dir';
% BC_dirichletLeft = BC;
% Lagrange multipliers to impose pressure jump with generalized periodicity
BC.dirichletLeft        = false;
BC.dirichletTopBot      = false;
BC.slidingTopBot        = false;
BC.periodic.topBot      = true;
BC.periodic.leftRight   = true;
BC.periodic.p.topBot    = true;
BC.periodic.p.leftRight = true;
BC.periodic.p.lagrange  = true;
BC.label = 'lagrangeJump';
BC_lagrangeJump = BC;
% Pressure gradient with periodicity on left-right pressure
BC.dirichletLeft        = false;
BC.dirichletTopBot      = false;
BC.slidingTopBot        = false;
BC.periodic.topBot      = true;
BC.periodic.leftRight   = true;
BC.periodic.p.topBot    = true;
BC.periodic.p.leftRight = true;
BC.periodic.p.lagrange  = false;
BC.label = 'pressGrad';
BC_pgrad = BC;
% Pressure gradient without periodicity on left-right pressure
BC.dirichletLeft        = false;
BC.dirichletTopBot      = false;
BC.slidingTopBot        = false;
BC.periodic.topBot      = true;
BC.periodic.leftRight   = true;
BC.periodic.p.topBot    = true;
BC.periodic.p.leftRight = false;
BC.periodic.p.lagrange  = false;
BC.label = 'pressGradNoPerX';
BC_pgrad_noperx = BC;