clc; close all; clear all;

%% --- Add paths ---
addpath('./input') 

%% --- Parameters ---
mu = 1e-3;          % [Pa.s] dynamic viscosity of fluid
alpha = 0.931;      % Drummond–Tahir model constant (example)
beta  = -1.476;     % Drummond–Tahir model constant (example)
pixel_size = 50e-9; % [m] per pixel (update from scale bar in image)
L = 1e-4;           % [m] characteristic length for Darcy number
C_kc = 5;           % Kozeny constant

%% --- Load and preprocess image ---
img = imread('mucus1.png');
if size(img,3) == 3
    img_gray = rgb2gray(img); % convert to grayscale if RGB
else
    img_gray = img;
end

% --- Segment pores and solids using Otsu threshold ---
level = graythresh(img_gray);  % automatic threshold
solidMask = imbinarize(img_gray, level); % 1 = solid fibers, 0 = pore space

% --- Compute packing fraction φ and porosity ε ---
phi = sum(solidMask(:)) / numel(solidMask); % solid fraction
epsilon = 1 - phi; % porosity

% --- Estimate characteristic radius R' from equivalent diameters ---
stats = regionprops(solidMask, 'EquivDiameter');
meanDiameter_pixels = mean([stats.EquivDiameter]); % mean object diameter
Rprime = (meanDiameter_pixels * pixel_size) / 2; % [m]

% --- Drummond–Tahir model ---
kprime_DT = (Rprime^2 / mu) * (1/(4*phi)) * ...
            (log(phi) + alpha + beta * phi^2 + 2*phi);
Da_DT = kprime_DT / L^2;

% --- Kozeny–Carman model ---
Sv = 2 / Rprime; % specific surface area [m^-1]
kprime_KC = (epsilon^3) / (C_kc * Sv^2 * (1 - epsilon)^2);
Da_KC = kprime_KC / L^2;

%% --- Display results ---
fprintf('--- Geometry parameters ---\n');
fprintf('Solid fraction φ: %.4f\n', phi);
fprintf('Porosity ε: %.4f\n', epsilon);
fprintf('Characteristic radius R'': %.4e m\n', Rprime);
fprintf('Specific surface area S_v: %.4e m^-1\n', Sv);

fprintf('\n--- Permeability results ---\n');
fprintf('DT model:  k'' = %.4e m^2, Da = %.4e\n', kprime_DT, Da_DT);
fprintf('KC model:  k'' = %.4e m^2, Da = %.4e\n', kprime_KC, Da_KC);

%% --- Visualise segmentation ---
figure;
subplot(1,2,1); imshow(img_gray); title('Original SEM image');
subplot(1,2,2); imshow(solidMask); title('Segmented solid phase');
