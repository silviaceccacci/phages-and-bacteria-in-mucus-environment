% =========================================================================
% 3D MUCIN SURFACE MESH GENERATOR (CONTINUOUS SWEPT TUBE)
% Generates watertight .stl and verifies the 3D actual density gradient
% =========================================================================
clc; clear all; close all;

addpath('./output/paraview')
%% --- 1. Domain Dimensions (in Nanometers) ---
M = 1000;       % X dimension
N = 1000;       % Y dimension
L_z = 15000;    % Z dimension (Height)

%% --- 2. Polymer Characteristics ---
strand_len = 500;   % Total length of one strand (nm)
kuhn_len = 30;      % Step size (nm)
radius_nm = 5.0;    % Cylinder radius (10nm thickness)
n_sides = 6;        % Number of sides on the cylinder
num_steps = floor(strand_len / kuhn_len); 

%% --- 3. Density Distribution & Strand Calculation ---
peak_mucin_fraction = 0.05; % Target 5% mucin density at the bottom
lambda = 2000;              % Decay constant (2 microns)

% 1. Mathematical Volume of a single continuous cylindrical strand
strand_vol = strand_len * (pi * radius_nm^2); 
% 2. Total required mucin volume to fill the exponential decay curve
target_vol = peak_mucin_fraction * (M * N * lambda);
% 3. Calculate total strands (with a 1.25 correction factor for overlaps)
correction_factor = 1.25; 
total_strands = round((target_vol / strand_vol) * correction_factor);

fprintf('Target Peak Density: %.1f%%\n', peak_mucin_fraction * 100);
fprintf('Generating exactly %d strands to reach target volume...\n', total_strands);

%% --- 4. Pre-allocate Arrays ---
% For the 3D STL Mesh
V_list = cell(total_strands * 5, 1); 
F_list = cell(total_strands * 5, 1);
shape_idx = 0;
v_offset = 0; 

% For the Density Verification Slicer (Memory Tracker)
max_total_steps = total_strands * num_steps;
segments = zeros(max_total_steps, 6);
segment_count = 0;

%% --- 5. Generate Polymer Mesh ---
for f = 1:total_strands
    % Birth Coordinates
    curr_x = rand() * M;
    curr_y = rand() * N;
    curr_z = -lambda * log(1 - rand() * (1 - exp(-L_z/lambda)));
    curr_z = max(0, min(L_z, curr_z)); 
    
    current_segment_pts = [curr_x, curr_y, curr_z];
    
    prev_u = [0, 0, 0]; 
    stiffness = 0.5; 
    
    for s = 1:num_steps
        % Random Walk Physics
        valid = false;
        while ~valid
            u = randn(1, 3); 
            u = u / norm(u); 
            if s == 1 || dot(u, prev_u) > stiffness
                valid = true; 
            end
        end
        prev_u = u;
        
        next_x = curr_x + kuhn_len * u(1);
        next_y = curr_y + kuhn_len * u(2);
        next_z = curr_z + kuhn_len * u(3);
        
        % ---> TRACK THE SEGMENT FOR DENSITY CALCULATION <---
        segment_count = segment_count + 1;
        segments(segment_count, :) = [curr_x, curr_y, curr_z, next_x, next_y, next_z];
        
        % Z-Boundary: Bounce (Reflect)
        if next_z < 0
            next_z = abs(next_z); u(3) = -u(3); prev_u = u;
        elseif next_z > L_z
            next_z = L_z - (next_z - L_z); u(3) = -u(3); prev_u = u;
        end
        
        % X/Y-Boundary: Periodic Check
        crossed_boundary = (next_x < 0 || next_x >= M || next_y < 0 || next_y >= N);
        
        if crossed_boundary
            % 1. Mesh the segment BEFORE it crosses the boundary
            if size(current_segment_pts, 1) > 1
                [V_tube, F_tube] = generate_swept_tube(current_segment_pts, radius_nm, n_sides);
                if ~isempty(V_tube)
                    shape_idx = shape_idx + 1;
                    V_list{shape_idx} = V_tube;
                    F_list{shape_idx} = F_tube + v_offset;
                    v_offset = v_offset + size(V_tube, 1);
                end
            end
            
            % 2. Teleport (Wrap) to the other side
            curr_x = mod(next_x, M);
            curr_y = mod(next_y, N);
            curr_z = next_z;
            
            % 3. Start a brand new segment
            current_segment_pts = [curr_x, curr_y, curr_z];
        else
            % Did not cross boundary: Keep appending to current segment
            curr_x = next_x; 
            curr_y = next_y; 
            curr_z = next_z;
            current_segment_pts(end+1, :) = [curr_x, curr_y, curr_z];
        end
    end
    
    % End of strand: Mesh the final remaining segment
    if size(current_segment_pts, 1) > 1
        [V_tube, F_tube] = generate_swept_tube(current_segment_pts, radius_nm, n_sides);
        if ~isempty(V_tube)
            shape_idx = shape_idx + 1;
            V_list{shape_idx} = V_tube;
            F_list{shape_idx} = F_tube + v_offset;
            v_offset = v_offset + size(V_tube, 1);
        end
    end
end

% Trim empty memory from segments
segments = segments(1:segment_count, :);

%% --- 6. Combine and Export as STL File ---
fprintf('Compiling watertight geometry data...\n');
All_Vertices = cell2mat(V_list(1:shape_idx));
All_Faces = cell2mat(F_list(1:shape_idx));

export = fullfile('./output/paraview', 'mucus_mesh_3D.stl');
TR = triangulation(All_Faces, All_Vertices);
stlwrite(TR, export);

fprintf('Success! Exported solid surface mesh with %d triangles', size(All_Faces, 1));

%% HELPER FUNCTION: GENERATE CONTINUOUS SWEPT TUBE (WITH END CAPS)
% =========================================================================
function [Vertices, Faces] = generate_swept_tube(points, radius, n_sides)
    num_pts = size(points, 1);
    if num_pts < 2
        Vertices = []; Faces = []; return;
    end
    
    T = zeros(num_pts, 3);
    for i = 1:num_pts
        if i == 1
            T(i,:) = points(2,:) - points(1,:);
        elseif i == num_pts
            T(i,:) = points(end,:) - points(end-1,:);
        else
            T(i,:) = points(i+1,:) - points(i-1,:); 
        end
        T(i,:) = T(i,:) / norm(T(i,:));
    end
    
    N_vec = zeros(num_pts, 3);
    B_vec = zeros(num_pts, 3);
    
    temp = cross(T(1,:), [0 0 1]);
    if norm(temp) < 1e-4, temp = cross(T(1,:), [0 1 0]); end
    N_vec(1,:) = temp / norm(temp);
    B_vec(1,:) = cross(T(1,:), N_vec(1,:));
    
    for i = 2:num_pts
        axis = cross(T(i-1,:), T(i,:));
        if norm(axis) > 1e-6
            axis = axis / norm(axis);
            angle = acos(max(-1, min(1, dot(T(i-1,:), T(i,:)))));
            K = [0 -axis(3) axis(2); axis(3) 0 -axis(1); -axis(2) axis(1) 0];
            R = eye(3) + sin(angle)*K + (1-cos(angle))*(K^2);
            N_vec(i,:) = (R * N_vec(i-1,:)')';
        else
            N_vec(i,:) = N_vec(i-1,:);
        end
        B_vec(i,:) = cross(T(i,:), N_vec(i,:));
    end
    
    theta = linspace(0, 2*pi, n_sides+1); theta(end) = [];
    Vertices = zeros(num_pts * n_sides + 2, 3); 
    
    v_idx = 1;
    for i = 1:num_pts
        for j = 1:n_sides
            Vertices(v_idx, :) = points(i,:) + radius * (cos(theta(j))*N_vec(i,:) + sin(theta(j))*B_vec(i,:));
            v_idx = v_idx + 1;
        end
    end
    
    cap1_idx = v_idx;
    cap2_idx = v_idx + 1;
    Vertices(cap1_idx, :) = points(1,:);
    Vertices(cap2_idx, :) = points(end,:);
    
    Faces = zeros((num_pts-1) * n_sides * 2 + n_sides * 2, 3); 
    f_idx = 1;
    for i = 1:(num_pts-1)
        ring1_start = (i-1) * n_sides;
        ring2_start = i * n_sides;
        
        for j = 1:n_sides
            next_j = mod(j, n_sides) + 1;
            p1 = ring1_start + j;
            p2 = ring1_start + next_j;
            p3 = ring2_start + j;
            p4 = ring2_start + next_j;
            
            Faces(f_idx, :) = [p1, p2, p3];
            Faces(f_idx+1, :) = [p2, p4, p3];
            f_idx = f_idx + 2;
        end
    end
    
    ring1_start = 0;
    ring2_start = (num_pts-1) * n_sides;
    for j = 1:n_sides
        next_j = mod(j, n_sides) + 1;
        Faces(f_idx, :) = [cap1_idx, ring1_start + next_j, ring1_start + j];
        f_idx = f_idx + 1;
        Faces(f_idx, :) = [cap2_idx, ring2_start + j, ring2_start + next_j];
        f_idx = f_idx + 1;
    end
end