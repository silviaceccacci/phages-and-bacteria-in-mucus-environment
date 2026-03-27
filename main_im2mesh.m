clc; clear all; close all;

addpath('./input')
addpath('./image2mesh/src')
addpath_tools()

%% Global variable to find adaptation tool (BAMG)
global machine;
%machine = 'abelm2';
%machine = 'silvia';
machine = 'jose';
do_adapt_mesh = false;

%% Input files
fileName = 'mucus1';
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
if(channels>1) 
   disp('Get first channel, b&w')
   Z = Z(:,:,1);
end

ximage = 0:(n-1);
yimage = (m-1):-1:0;
[x,y] = meshgrid(ximage,yimage);
    
X = [x(:), y(:)];
[T] = delaunay(X(:,1),X(:,2));

Z = Z(:);

iexport = iexport+1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = Z;
if(do_exportInitialImage) 
    exportMeshParaview(X,T,options)
end
%% Adapted mesh
threshold_mucus = mean(Z)*percentage_image_to_consider_mucin;
fact_tanh = 100;
Z = Z-threshold_mucus;
Z = tanh(Z*fact_tanh);

total_area = sum(ones(size(Z)));
mucin_area = sum((1+Z)/2);
density_mucin = mucin_area/total_area

iexport = iexport+1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = Z ;
exportMeshParaview(X,T,options)

if(do_adapt_mesh) 
    adapt_variable.f = 1;
    sol.f = Z;
    [Xh,Th,sol] = adapt_multipleSol(X,T,data,adapt_variable,sol,doPeriodic);
    Zh = sol.f;
else
    Xh=X;
    Th=T;
    Zh=Z;
end

iexport = iexport+1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = Zh;
exportMeshParaview(Xh,Th,options)

fprintf('Num nodes: %d \n',size(Xh,1))
fprintf('Num elems: %d \n',size(Th,1))
%% Force periodicity
if(doPeriodic)
    mesh.X = Xh;
    mesh.T = Th;
    mesh= makePeriodicMesh(mesh,hmin);
    isPeriodic = checkPeriodicity(mesh);
end
%% Set reasonable permeability
% Interpolate using bilinear method
Zh = interp2(x, y, img, mesh.X(:,1), mesh.X(:,2), 'linear');
% Transform:
Zh = Zh-threshold_mucus;
Zh = tanh(Zh*fact_tanh);

perme = Zh; % here in range -1,1
perme = (1-perme)/2; % now zero for mucin, 1 rest

iexport = iexport+1;
options.exportName = [fileName '_' int2str(iexport)];
options.f = Zh;
exportMeshParaview(mesh.X,mesh.T,options)
%% Save mesh
x_max = max(mesh.X(:,1));
mesh.X = mesh.X*scale/x_max;

mesh.perme = perme;
save(['./output/' fileName '_mesh'],"mesh")
%text_density = int2str(round(density_mucin*100));
%save(['./out/' fileName '_mesh'],"mesh")

