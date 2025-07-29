function [domain,mesh,model,parameters]=adimensionalize(domain,mesh,model,parameters)

    % DOMAIN
%     x0 = [domain.x_LL domain.y_LD];
%     x1 = [domain.x_LR domain.y_LU];
%     domain.Lchar= norm(x1-x0);
    domain.Lchar=1e-6; % 1 micrometer, to make it constant throught different images?
    domain.origin = [domain.x_LL,domain.y_LD];  
    
    mesh.X(:,1) = mesh.X(:,1)-domain.origin(1);
    mesh.X(:,2) = mesh.X(:,2)-domain.origin(2);
    mesh.X = mesh.X/domain.Lchar;
    
    domain.x_LL = min(mesh.X(:,1));
    domain.x_LR = max(mesh.X(:,1));
    domain.y_LD = min(mesh.X(:,2));
    domain.y_LU = max(mesh.X(:,2));
    domain.Lx = domain.x_LR-domain.x_LL;
    domain.Ly = domain.y_LU-domain.y_LD;
    mesh.domain = domain;
    fprintf('Domain adim: [%4.2e %4.2e]x[%4.2e %4.2e]\n',domain.x_LL,domain.x_LR,domain.y_LD,domain.y_LU)
    
    disp('ENSURE THAT PRESSURE GRAD/JUMP IS CORRECT!!!!!!')
    % adim pressure grad:
    %parameters.pressureGrad = parameters.pressureGrad/domain.Lchar;
    Pref = domain.Lchar/( model.viscosity*parameters.umag_in );
    parameters.pressureGrad = parameters.pressureGrad/Pref; 
    % compute total jump:
    parameters.pressureJump = parameters.pressureGrad*domain.Lx;    
    
    % Boundary conditions
    parameters.umag_in_dimensional = parameters.umag_in;
    parameters.umag_in = 1.0;

    % PERME AND VISCOSITY
    model.viscosity_dim  = model.viscosity ;
    model.viscosity  = 1.0 ;
    %model.perme = model.perme / (domain.Lchar)^2;
    % lets consider a non dimensional input perme (easier)
    % this adim perme is the Darcy number

    model.u_dir      = @(X)(bc_const(X,parameters));
    %model.source     = @(X)(source_zero(X,parameters));
    model.source = @(X)([-parameters.pressureGrad; 0.0]);
