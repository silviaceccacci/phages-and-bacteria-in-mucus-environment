clc; clear all; close all;

addpath('./input')

% Modular Mucus Surface Generator
% Methodology: Pérez-Ràfols & Almqvist (2019)
clear; clc;

%% --- 0. SELECTION SETTINGS ---
% Choose PSD: 1 = Fixed Ring, 2 = Aberrant Multi-scale, 3 = Fractal
% (Self-Affine), 4 = Fixed Ring with Noise
psd_type = 2;

% Choose HPD: 1 = Weibull, 2 = Bimodal (Min-Gaussian), 3 = Gaussian
hpd_type = 3;
%% 1. Grid & Physical Scaling
pixel_size = 5;
width_nm = 5500;
height_nm = 4300;

M = round(width_nm / pixel_size);  % Number of points in x
N = round(height_nm / pixel_size); % Number of points in y

target_pore_nm = 180;   % Primary target spacing

%% 2. Define Target Power Spectrum (PSD)
% Use M and N separately for the meshgrid to create a rectangular spectral domain
[kx, ky] = meshgrid((-M/2:M/2-1), (-N/2:N/2-1));

% Calculate frequencies based on specific dimensions 
% Note: k is now an NxM matrix
k = sqrt((kx * (2*pi/width_nm)).^2 + (ky * (2*pi/height_nm)).^2);
k_peak = 2*pi / target_pore_nm;
sigma_k = k_peak/2;

switch psd_type
    case 1 % Standard 180nm Ring
        PSD_target = exp(-(k - k_peak).^2 / (2 * (sigma_k)^2));
    case 2 % Aberrant Multi-scale
        k_clump = 2*pi / 500; 
        sigma_c = k_clump/2;
        PSD_target = 0.5 * exp(-(k - k_peak).^2 / (2 * (sigma_k)^2)) + ... 
                     1.5 * exp(-(k - k_clump).^2 / (2 * (sigma_c)^2));
    case 3 % Fractal / Self-Affine
        H_f = 0.8; 
        PSD_target = k.^(-2*(1 + H_f));

    case 4 % Stochastic 180nm Ring (Jagged Frontiers)
        A = 1.0;            % Strength of the ring structure
        noise_floor = 0.008; % The "Aberration" factory
        
        % Create the standard Ring
        PSD_ring = A * exp(-(k - k_peak).^2 / (2 * sigma_k^2));
        
        % Add the White Noise Floor
        PSD_target = PSD_ring + noise_floor;
        
end

% Set DC component to zero to ensure zero mean
PSD_target(N/2+1, M/2+1) = 0; 

%% 3. Define Target HPD
% Ensure we generate exactly N*M points

switch hpd_type
    case 1 % Weibull HPD
        lambda = 1; % keep it at 1 for standard visualization
        k = 1.2; % k = 3.6 nearly gaussian
        target_hpd_values = wblrnd(lambda, k, [N*M, 1]);
        
    case 2 % Bimodal HPD
        % PARAMETERS TO CHANGE:
        c_val = 5;          % Distance between peaks (Try 5 vs 15)
        mucin_fraction = 0.7; % Mucin fraction
        
        % Create two separate populations
        num_mucin = round(N*M * mucin_fraction);
        num_fluid = (N*M) - num_mucin;
        
        % Peak 1 (Fluid) centered at 0
        fluid_peak = randn(num_fluid, 1);
        
        % Peak 2 (Mucin) centered at c_val
        mucin_peak = randn(num_mucin, 1) + c_val;
        
        % Combine and sort
        target_hpd_values = sort([fluid_peak; mucin_peak]);
        
    case 3 % Gaussian HPD (Standard/Smooth)
        mu = 0;      % Mean
        sigma = 1.5; % Standard Deviation
        
        % This generates N*M values following N(mu, sigma)
        target_hpd_values = mu + sigma .* randn(N*M, 1);
end

% Sort the values for the rank-ordering algorithm
target_hpd_values = sort(target_hpd_values);

%% 4. Iterative Algorithm
% Initialize with a seed matching the PSD
phase = 2 * pi * rand(N, M);
Z_f = sqrt(PSD_target) .* exp(1i * phase);
z_current = real(ifft2(ifftshift(Z_f)));

for iter = 1:15
    % Step A: Force HPD via Rank Ordering
    [~, idx] = sort(z_current(:));
    z_new = zeros(N*M, 1);
    z_new(idx) = target_hpd_values; 
    z_hpd = reshape(z_new, [N, M]);
    
    % Step B: Force PSD via Spectral Correction
    Z_f_current = fftshift(fft2(z_hpd));
    Z_f_corrected = sqrt(PSD_target) .* exp(1i * angle(Z_f_current));
    z_current = real(ifft2(ifftshift(Z_f_corrected)));
end

%% 5. Final Scaling & Export
mucin_image = (z_hpd - min(z_hpd(:))) / (max(z_hpd(:)) - min(z_hpd(:)));
mucin_image = uint8(255 * mucin_image);

%% 6. Visualization
figure('Position', [100, 100, 1100, 400]);
subplot(1,2,1); imshow(mucin_image);
title(sprintf('Mucin Mesh (PSD: %d, HPD: %d)', psd_type, hpd_type));
subplot(1,2,2); histogram(mucin_image, 50, 'FaceColor', 'k');
title('Obtained Height Probability Distribution (HPD)');

% Save the mucin topography image to a file
routeimg = fullfile('input', 'PSD_topography.png');
imwrite(mucin_image, routeimg);