clc; clear all; close all;

addpath('./input')
addpath('./image2mesh/src')
addpath_tools()

%% Global variable to find adaptation tool (BAMG)
global machine;
%machine = 'abelm2';
%machine = 'silvia';
machine = 'jose';
do_adapt_mesh = true;

%% Input files
fileName = 'PSD_topography';
format = 'png';

doPeriodic = true;

doPartOfImage = false;
portion = 0.25;

micron = 1e-6; % Conversion for micrometers
nano = 1e-9; % Conversion for nanometer
x_min = 0 * micron;
y_min = 0 * micron;
x_max  = 5.5 * micron; 
y_max = 3.5 * micron; 
scale = x_max;

if(doPeriodic)
    periodic_percent = 0.005;%0.005
    scale = scale*(1+periodic_percent);
end

do_exportInitialImage = false;
%% Physical parameters
disp('Threshold to consider what is mucin and not is relevant...')
percentage_image_to_consider_mucin = 0.9;

%% Mesh parameters
[data] = parameters_bamg();

hmin = 1.0;
hmax = 8.0;
data.hmin   = num2str(hmin) ;%'2.0' ;
data.hmax   = num2str(hmax) ;%'10.0' ;
data.anisomax = '2';     %0 without anisomax
data.ratio = '1.1';
%% Read image
iexport = 0;

imageName = [fileName '.' format];

img = imread(imageName); 
img = im2double(img); % Normalize the image values between 0 and 1

% RGB = double(imread(imageName)); 
img = img(:,:,1);

img = img(1:end-5,5:1050,:); % to remove croped things

fprintf('Size of image: %d %d\n',size(img))

if(doPartOfImage) 
    n = min([size(img,1),size(img,2)]);
    n = ceil(n*portion);
    img = img(1:n,1:n,:);
    %I = I(1:50,1:n,:);
end

figure; imshow(img); title('Original')
%% Periodic image
if(doPeriodic)
    %img = [img  img(:,end:-1:1)];
    %img = [img;  img(end:-1:1,:)];
    
    nx = size(img,2);
    nx_per = ceil(nx*periodic_percent);
    nx_per = max(nx_per,5);

    im_out = img(:, end);
    im_in = img(:, 1);
    lambda = (1:nx_per)/nx_per;
    im_blend = im_out*(1-lambda) + im_in*lambda;
    img = [img  im_blend];
    

    im_bot = img(end, :);
    im_top = img(1, :);
    lambda = (1:nx_per)/nx_per;
    lambda = lambda';
    im_blend = (1-lambda)*im_bot + lambda*im_top;
    img = [img ; im_blend];
    
%     img = [img  img];
%     img = [img ; img];
    figure; imshow(img); title('Periodic')
end
%% Filter
%target_density = 0.1; % Desired percentage (0-4%)
%[img]=filter_for_density(img,target_density);

% total_area = sum(ones(size(Z)));
% threshold_mucus = mean(Z);%/10;
% f = Z-threshold_mucus;
% f = tanh(f*fact_tanh);
% f = (f+1)/2;
% % while(is_desired_density == false)
% %     threshold_mucus = threshold_mucus*1.1;
% %     f = Z-threshold_mucus;
% %     f = tanh(f*fact_tanh);
% %     f = (f+1)/2;
% %     mucin_area = sum(f);
% %     density_mucin = mucin_area/total_area;
% % 
% %     is_desired_density = abs(target_density-density_mucin)<1e-2;
% % end
% % density_mucin
%% Structured mesh
Z = img;
[m, n, channels] = size(Z);
if channels > 1
    disp('Get first channel, b&w')
    Z = Z(:,:,1);
end

%% Thresholding and tanh smoothing
threshold_mucus = max(Z(:)) * 0.65; % threshold to consider what is mucin and what is not
fact_tanh = 100;
Z_smooth = tanh((double(Z) - threshold_mucus) * fact_tanh);  % 2D, values in [-1, 1]

total_area  = numel(Z_smooth);
mucin_area  = sum((1 + Z_smooth(:)) / 2);  % mucin ~ +1 side
density_mucin = mucin_area / total_area;
fprintf('Density of the mucin: %.4f\n', density_mucin);

is_fluid = Z_smooth < 0;   % 2D logical mask

%% Build pixel grid
ximage = 0:(n-1);
yimage = (m-1):-1:0;
[x, y] = meshgrid(ximage, yimage);

X_all      = [x(:), y(:)];    % all pixel coordinates
fluid_mask = is_fluid(:);      % logical column vector

%% Export initial full image for reference
if do_exportInitialImage
    T_all = delaunay(X_all(:,1), X_all(:,2));
    iexport = iexport + 1;
    options.exportName = [fileName '_' int2str(iexport)];
    options.f = Z_smooth(:);           % full grid, matches X_all
    exportMeshParaview(X_all, T_all, options)
end

%% Restrict to fluid nodes only
X = X_all(fluid_mask, :);
T = delaunay(X(:,1), X(:,2));

%% Remove triangles bridging across mucin (centroid check)
centroids = (X(T(:,1),:) + X(T(:,2),:) + X(T(:,3),:)) / 3;

centroid_fluid = interp2(x, y, double(is_fluid), ...
                         centroids(:,1), centroids(:,2), 'nearest', 0);

T = T(centroid_fluid > 0.5, :);   % keep only fluid triangles

used_nodes = unique(T(:));
node_map   = zeros(size(X, 1), 1);
node_map(used_nodes) = 1:numel(used_nodes);

X = X(used_nodes, :);
T = node_map(T);

fluid_indices        = find(fluid_mask);
fluid_mask           = false(size(fluid_mask));
fluid_mask(fluid_indices(used_nodes)) = true;
%% Export fluid-only mesh (pre-adaptation)
iexport = iexport + 1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = Z_smooth(fluid_mask);      % fluid nodes only, matches X
exportMeshParaview(X, T, options)

%% Identify boundary edges
all_edges  = [T(:,[1,2]); T(:,[2,3]); T(:,[1,3])];
all_edges  = sort(all_edges, 2);
[uniq_edges, ~, ic] = unique(all_edges, 'rows');
edge_count = accumarray(ic, 1);
boundary_edges = uniq_edges(edge_count == 1, :);

on_image_border = ( X(boundary_edges(:,1), 1) == 0     | ...
                    X(boundary_edges(:,1), 1) == (n-1) | ...
                    X(boundary_edges(:,1), 2) == 0     | ...
                    X(boundary_edges(:,1), 2) == (m-1) | ...
                    X(boundary_edges(:,2), 1) == 0     | ...
                    X(boundary_edges(:,2), 1) == (n-1) | ...
                    X(boundary_edges(:,2), 2) == 0     | ...
                    X(boundary_edges(:,2), 2) == (m-1) );

solid_wall_edges   = boundary_edges(~on_image_border, :);
domain_inlet_edges = boundary_edges( on_image_border, :);

%% Adapted mesh refinement driven by distance to mucin wall
dist_to_mucin = bwdist(~is_fluid);
dist_fluid    = double(dist_to_mucin(fluid_mask));

% Define the narrow band where BAMG is allowed to see curvature
band_width = 10;  % in pixels

% Linear ramp clipped at band_width 
dist_clipped = min(dist_fluid, band_width);

% Cosine taper: smooth from 0 at wall to 1 at band edge, then flat
% Hessian is nonzero only very close to wall, zero at band edge and beyond
refine_field = 0.5 * (1 - cos(pi * dist_clipped / band_width));

adapt_variable.f = 1;
sol.f = refine_field;
[Xh, Th, sol] = adapt_multipleSol(X, T, data, adapt_variable, sol, doPeriodic);
Zh = sol.f;

%% Export adapted fluid mesh
iexport = iexport + 1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = Zh;
exportMeshParaview(Xh, Th, options)

%% Calculate new count of boundary edges and compare it to pre-adaptation
all_edges_h    = [Th(:,[1,2]); Th(:,[2,3]); Th(:,[1,3])];
all_edges_h    = sort(all_edges_h, 2);
[uniq_edges_h, ~, ic_h] = unique(all_edges_h, 'rows');
edge_count_h   = accumarray(ic_h, 1);
boundary_edges_h = uniq_edges_h(edge_count_h == 1, :);

on_image_border_h = ( Xh(boundary_edges_h(:,1), 1) == 0     | ...
                      Xh(boundary_edges_h(:,1), 1) == (n-1) | ...
                      Xh(boundary_edges_h(:,1), 2) == 0     | ...
                      Xh(boundary_edges_h(:,1), 2) == (m-1) | ...
                      Xh(boundary_edges_h(:,2), 1) == 0     | ...
                      Xh(boundary_edges_h(:,2), 1) == (n-1) | ...
                      Xh(boundary_edges_h(:,2), 2) == 0     | ...
                      Xh(boundary_edges_h(:,2), 2) == (m-1) );

solid_wall_edges_h   = boundary_edges_h(~on_image_border_h, :);
domain_inlet_edges_h = boundary_edges_h( on_image_border_h, :);

fprintf('--- Pre-adaptation ---\n')
fprintf('  Boundary edges (domain):       %d\n', size(domain_inlet_edges, 1))
fprintf('  Solid wall edges (mucin wall): %d\n', size(solid_wall_edges, 1))
fprintf('--- Post-adaptation ---\n')
fprintf('  Boundary edges (domain):       %d\n', size(domain_inlet_edges_h, 1))
fprintf('  Solid wall edges (mucin wall): %d\n', size(solid_wall_edges_h, 1))

fprintf('Num nodes: %d\n', size(Xh, 1))
fprintf('Num elems: %d\n', size(Th, 1))
fprintf('Density of mucin (solid fraction): %.4f\n', density_mucin)

%% Force periodicity
if(doPeriodic)
    mesh.X = Xh;
    mesh.T = Th;
    mesh = makePeriodicMeshStokes(mesh,hmin);
    isPeriodic = checkPeriodicity(mesh);
end

%% Export per periodic fluid mesh
iexport = iexport + 1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = ones(size(mesh.X, 1), 1); % just to check if the geometry is correct
exportMeshParaview(mesh.X, mesh.T, options)

%% Save mesh
x_max = max(mesh.X(:,1));
mesh.X = mesh.X*scale/x_max;

save(['./output/' fileName '_mesh'],"mesh")
%text_density = int2str(round(density_mucin*100));
%save(['./out/' fileName '_mesh'],"mesh")

