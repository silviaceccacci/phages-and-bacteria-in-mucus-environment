% =========================================================================
% 3D SOLID MUCIN VOLUME GENERATOR, DENSITY CHECK & 2D SLICER
% Generates a solid 3D logical matrix (1 = Mucin, 0 = Fluid)
% =========================================================================
clc; clear all; close all;

addpath('./input')

%% --- 1. Domain Dimensions & Voxel Resolution ---
M = 1000;       % X dimension (nm)
N = 1000;       % Y dimension (nm)
L_z = 15000;    % Z dimension (nm)

voxel_size = 2.5; % 1 voxel = 1 nm (Adjust this based on your mesher's needs)
grid_X = round(M / voxel_size);
grid_Y = round(N / voxel_size);
grid_Z = round(L_z / voxel_size);

%% --- 2. Polymer Characteristics & Density ---
strand_len = 500;   
kuhn_len = 30;      
radius_nm = 5.0;    % Cylinder radius (10 nm total thickness)
r_vox = radius_nm / voxel_size; % Radius in voxel units

num_steps = floor(strand_len / kuhn_len); 
peak_mucin_fraction = 0.05; 
lambda = 2000;              

% Calculate strands (Using a correction factor to overcome bleed/overlap)
strand_vol = strand_len * (pi * radius_nm^2); 
target_vol = peak_mucin_fraction * (M * N * lambda);
correction_factor = 1; % Increased to ensure actual 3D volume hits 5%
total_strands = round((target_vol / strand_vol) * correction_factor);

fprintf('Generating %d strands to build the solid volume...\n', total_strands);

%% --- 3. Run the Random Walk & Track Segments ---
max_total_steps = total_strands * num_steps;
segments = zeros(max_total_steps, 6);
segment_count = 0;

for f = 1:total_strands
    curr_x = rand() * M;
    curr_y = rand() * N;
    % Using exponential distribution for density gradient along Z
    curr_z = -lambda * log(1 - rand() * (1 - exp(-L_z/lambda)));
    curr_z = max(0, min(L_z, curr_z)); 
    
    prev_u = [0, 0, 0]; 
    stiffness = 0.5; 
    
    for s = 1:num_steps
        valid = false;
        while ~valid
            u = randn(1, 3); u = u / norm(u); 
            if s == 1 || dot(u, prev_u) > stiffness, valid = true; end
        end
        prev_u = u;
        
        next_x = curr_x + kuhn_len * u(1);
        next_y = curr_y + kuhn_len * u(2);
        next_z = curr_z + kuhn_len * u(3);
        
        segment_count = segment_count + 1;
        segments(segment_count, :) = [curr_x, curr_y, curr_z, next_x, next_y, next_z];
        
        % Bounce Z
        if next_z < 0, next_z = abs(next_z); u(3) = -u(3); prev_u = u;
        elseif next_z > L_z, next_z = L_z - (next_z - L_z); u(3) = -u(3); prev_u = u; end
        
        % Wrap X and Y
        curr_x = mod(next_x, M);
        curr_y = mod(next_y, N);
        curr_z = next_z;
    end
end
segments = segments(1:segment_count, :);

%% --- 4. Rasterize Solid Cylinders into 3D Volume Matrix ---
fprintf('Burning solid cylinders into %dx%dx%d Voxel Grid...\n', grid_X, grid_Y, grid_Z);

% Initialize the totally empty fluid domain (0s)
Volume_3D = false(grid_Y, grid_X, grid_Z);

% Create a small local 3D brush (a solid sphere) to speed up drawing
brush_size = ceil(r_vox);
[bx, by, bz] = meshgrid(-brush_size:brush_size, -brush_size:brush_size, -brush_size:brush_size);
brush_mask = (bx.^2 + by.^2 + bz.^2) <= r_vox^2;

% Draw every segment into the 3D matrix
for i = 1:size(segments, 1)
    p1 = segments(i, 1:3) / voxel_size;
    p2 = segments(i, 4:6) / voxel_size;
    
    % Interpolate points along the segment (step size = half a voxel)
    dist = norm(p2 - p1);
    num_pts = max(2, ceil(dist * 2));
    
    interp_x = linspace(p1(1), p2(1), num_pts);
    interp_y = linspace(p1(2), p2(2), num_pts);
    interp_z = linspace(p1(3), p2(3), num_pts);
    
    for pt = 1:num_pts
        % Center of the brush
        cx = round(interp_x(pt));
        cy = round(interp_y(pt));
        cz = round(interp_z(pt));
        
        % Handle periodic wrapping for X and Y
        cx = mod(cx - 1, grid_X) + 1;
        cy = mod(cy - 1, grid_Y) + 1;
        
        % Check Z bounds (it doesn't wrap, just stops at the floor/ceiling)
        if cz < 1 || cz > grid_Z, continue; end
        
        % Paint the solid brush onto the 3D volume
        for ix = -brush_size:brush_size
            for iy = -brush_size:brush_size
                for iz = -brush_size:brush_size
                    if brush_mask(iy + brush_size + 1, ix + brush_size + 1, iz + brush_size + 1)
                        % Wrapped coordinates for the brush edges
                        px = mod(cx + ix - 1, grid_X) + 1;
                        py = mod(cy + iy - 1, grid_Y) + 1;
                        pz = cz + iz;
                        
                        if pz >= 1 && pz <= grid_Z
                            Volume_3D(py, px, pz) = true;
                        end
                    end
                end
            end
        end
    end
end
fprintf('Done! The matrix "Volume_3D" is generated.\n');

%% --- 5. Calculate and Verify Z-Axis Density Distribution ---
fprintf('Calculating actual mucin density along the Z-axis...\n');

% 1. Count the number of mucin voxels (trues) in every XY plane slice
mucin_voxels_per_slice = squeeze(sum(sum(Volume_3D, 1), 2)); 

% 2. Divide by the total number of voxels in an XY plane to get volume fraction
actual_density = mucin_voxels_per_slice / (grid_X * grid_Y);

% 3. Create a Z-axis array in nanometers for plotting
z_array_nm = (1:grid_Z)' * voxel_size;

% 4. Calculate the theoretical target density
%    (Peak density decaying exponentially with constant lambda)
theoretical_density = peak_mucin_fraction * exp(-z_array_nm / lambda);

% Plot the actual vs theoretical density
figure('Name', 'Mucin Density Profile', 'Color', 'w');
plot(z_array_nm, actual_density, 'b-', 'LineWidth', 2);
hold on;
plot(z_array_nm, theoretical_density, 'r--', 'LineWidth', 2);
grid on;
xlabel('Z - Height (nm)', 'FontWeight', 'bold');
ylabel('Volume Fraction (Density)', 'FontWeight', 'bold');
title('Mucin Density Distribution along Z-axis', 'FontWeight', 'bold');
legend('Actual Generated Volume', 'Theoretical Target (Exponential)', 'Location', 'northeast');

%% --- 6. 2D Slicing (XZ Plane) and Pure PNG Save ---
fprintf('Extracting XZ slice and saving pure image as PNG...\n');

% Define the target plane and the slice index
slice_plane = 'XZ'; 
slice_index = round(grid_Y / 2); 
slice_index = min(max(slice_index, 1), grid_Y); 

% Extract the XZ plane directly from the 3D volume
% In a (Y, X, Z) coordinate system, fixing the first dimension (Y) yields the XZ plane.
% The ' operator (transpose) is used to align the orientation correctly.
slice_2D = squeeze(Volume_3D(slice_index, :, :))'; 

% Flip the image vertically so that Z=1 (bottom of the volume) appears at the bottom of the image
slice_2D_img = flipud(slice_2D);

% Generate the file path and save the logical matrix as a PNG (0 = black, 1 = white)
slice_png_path = fullfile('input', 'slice_2D_XZ.png');
imwrite(slice_2D_img, slice_png_path);

% Optional: Show a clean preview window
fig_slice = figure('Name', '2D XZ Mucin Slice Preview', 'Color', 'k'); 
imshow(slice_2D_img);
title(sprintf('XZ Slice Preview (Y = %d)', slice_index), 'Color', 'w');