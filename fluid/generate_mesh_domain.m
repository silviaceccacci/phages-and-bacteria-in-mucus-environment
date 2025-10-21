function [mesh,model,domain]=generate_mesh_domain(do_mesh_from_image,model,parameters)
if(do_mesh_from_image)
    meshName = [parameters.case_name '_mesh'];
    load(['./output/' meshName])
    fprintf('Num nodes: %d\n',size(mesh.X,1))
    
    % Set domain limits:
    domain.x_LL = min(mesh.X(:,1));
    domain.x_LR = max(mesh.X(:,1));
    domain.y_LD = min(mesh.X(:,2));
    domain.y_LU = max(mesh.X(:,2));
    mesh.domain = domain;

    fprintf('Domain:      [%4.2e %4.2e]x[%4.2e %4.2e]\n',domain.x_LL,domain.x_LR,domain.y_LD,domain.y_LU)
 
    % Set permeability:
    perme = mesh.perme; % from image

    min_p = min(perme);
    max_p = max(perme);
    minDarcyNum=parameters.minDarcyNum;
    maxDarcyNum=parameters.maxDarcyNum;
    factPerme = (maxDarcyNum-minDarcyNum)/(max_p-min_p);
    perme = minDarcyNum + (perme-min_p)*factPerme;
    model.perme = perme;
    clear perme;
else
    parameters.h_ini = 0.02;
    domain.LLcorner = [-1;-1];
    domain.URcorner = [1 ; 1];
    [mesh] = generateMesh(domain,parameters);

    %model.perme      = @perme_analytical;
    model.perme      = perme_analytical(mesh.X);
end