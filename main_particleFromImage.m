clear all; close all; clc;
%% Add paths
addpath('./particles') 
addpath('./fluid')
addpath('./fluid/solverNS/')
addpath('./interpolation')
addpath('./output')
%% TODOs:
disp('--------------------------------------------------------------------------------------------------------------')
disp('TODO')
disp(' - Particles solver:')
disp('--------------------------------------------------------------------------------------------------------------')
%% Define the folder to save all the files
outputFolder = './output/';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end
%% Domain parameters
disp('----> Read domain parameters')
micron = 1e-6; 
nano = 1e-9; 

case_name = 'mucus1';
load([case_name '_mesh'])
x_min = min(mesh.X(:,1));
x_max = max(mesh.X(:,1));
y_min = min(mesh.X(:,2));
y_max = max(mesh.X(:,2));
clear mesh;

load([case_name '_uInterp_1e-5'])

Omega_X = x_max - x_min; 
Omega_Y = y_max - y_min; 
domain = [Omega_X, Omega_Y];
%% Numerical parameters
disp('----> Read numerical and physical parameters')
dt = 1e-4; 
num_steps = 100; %10000;
%% Physical parameters
mu_water = 10^(-3);     % Dynamic viscosity (Pa s)
rho_water = 10^3;       % Mass density of phage (kg/m^3)
kB = 1.38e-23;          % Boltzmann constant (J/K)
T = 300;                % Temperature (Kelvin)
%% Phages
disp('----> Initialise phages')
rP = 100 * nano;                % Radius of phage (m)
rhoP = 10^5*rho_water;          % Mass density of phage (kg/m^3)
d_enc1 = rP;                    % Encounter distance for pha-bac attachemnt (lyse)
d_enc2 = rP/2;                  % Encounter distance for pha-bacInCl attachemnt (lyse)
d_enc3 = rP/2;                  % Encounter distance for pha-COMcl attachemnt (lyse)
num_phages = 400;
for i = 1:num_phages
    phages(i) = Phage(rP, rhoP, mu_water, kB, T, domain, x_min, y_min, x_max, y_max);
    phages(i).id = i;
end
%% Bacteria
disp('----> Initialise bacteria')
lB = 3 * micron;                        % Length of bacterium (m)
wB = 0.5 * micron;                      % Width of bacterium (m)
rhoB = 5*10^3*rho_water;                % Mass density of bacterium (kg/m^3)
vB = 30 * micron;                       % Velocity of bacterium in run phase (m/s)
omega_T = 0.5;                          % Tumble frequency (1/s)
epsilon = 10 * kB * T;                  % Phage-bacteria iteraction strength
crit_distance_bacteria = 1 * micron;    % Encounter distance for bacteria-bacteria attachemnt (biofilm formation)
num_bacteria = 15;
max_num_bacteria = num_bacteria;
for i = 1:num_bacteria
    bacteria(i) = Bacterium(lB, wB, rhoB, mu_water, kB, T, dt, vB, omega_T, domain, ...
        x_min, y_min, x_max, y_max);
    bacteria(i).id = i;
    bacteria(i).phages_ids = [];
end
%% Clusters
disp('----> Initialise clusters at initial configuration')
D = 1e-12;
threshold = compute_clustering_threshold(U_interp, x_min, x_max, y_min, y_max, lB, wB, dt, D);
%threshold = 0.5 * micron;
[clusters, bacteria] = Cluster.form_clusters(bacteria, threshold, domain, x_min, y_min, x_max, y_max);

% Print formed clusters
fprintf('Formed %d clusters\n', length(clusters));
for i = 1:length(clusters)
    fprintf('Cluster %d: %d bacteria\n', i, clusters(i).size);
    ids = arrayfun(@(b) b.id, clusters(i).bacteria);
    fprintf('   Bacteria IDs: %s\n', mat2str(ids))
end

for i = 1:length(clusters) 
    clusters(i) = clusters(i).update_center_of_mass_and_group_velocity();
    clusters(i).id = i;
end
%% Inilitalise for trajectories
fig1 = figure(1); hold on;
positionsB_InCluster = []; 
comC = [];
for i = 1:length(clusters)
    for j = 1:length(clusters(i).bacteria)
        positionsB_InCluster(:, end+1) = clusters(i).bacteria(j).position(:);
    end
    if clusters(i).size > 1
        comC(:, end+1) = clusters(i).position;
    end
end

all_bacteria = [clusters.bacteria]; 
positionsB_InCluster = reshape([all_bacteria.position], 2, []);
xC = positionsB_InCluster(1, :); %x of bacteria in cluster
yC = positionsB_InCluster(2, :); %y of bacteria in cluster
h1 = plot(xC / micron, yC / micron, 'o', 'MarkerFaceColor', [0 0.4470 0.7410], ...
  'MarkerEdgeColor', [0 0 0], 'MarkerSize', 32);

all_clusters = [clusters.position];
comC = reshape(all_clusters, 2, []);
xCOM = comC(1, :); %x of com
yCOM = comC(2, :); %y of com
h2 = plot(xCOM / micron, yCOM / micron, 'x', 'MarkerFaceColor', [0 0 0], ...
  'MarkerEdgeColor', [0 0 0], 'MarkerSize', 32);

positionsB = reshape([bacteria.position], 2, []);
xB = positionsB(1, :); %x of bacteria not in cluster
yB = positionsB(2, :); %y of bacteria not in cluster
h3 = plot(xB / micron, yB / micron, 'o', 'MarkerFaceColor', [0.9290 0.6940 0.1250], ...
  'MarkerEdgeColor', [0 0 0], 'MarkerSize', 32);

positionsP = reshape([phages.position], 2, []);
xP = positionsP(1, :); %x of phages
yP = positionsP(2, :); %y of phages
h4 = plot(xP / micron, yP /micron, 'o', 'MarkerFaceColor', [0.8500 0.3250 0.0980], ...
  'MarkerEdgeColor', [0 0 0], 'MarkerSize', 8);

xlabel('$\Omega_x \ (\mu m)$', 'Interpreter', 'LaTeX', 'FontSize', 16);
ylabel('$\Omega_y \ (\mu m)$', 'Interpreter', 'LaTeX', 'FontSize', 16);
title('Phages, bacteria and clusters dynamics', 'Interpreter', 'LaTeX', 'FontSize', 16);
xlim([0 Omega_X / micron]);
ylim([0 Omega_Y / micron]);
set(gca, 'FontSize', 16);
time_text1 = text(0.05 * Omega_X / micron, 0.95 * Omega_Y / micron, 'Time: 0.0 s', 'FontSize', 12, 'Color', 'k');

videoFile1 = fullfile('output', 'particles_dynamics.avi');
video1 = VideoWriter(videoFile1);
video1.FrameRate = 10;
open(video1);
%% Allocation of variables
disp('----> Allocate variables')
coordP_over_time = zeros(num_steps, num_phages*2);
coordB_over_time = zeros(num_steps, max_num_bacteria*2);
max_clusters = length(bacteria); 
coordC_over_time = zeros(num_steps, max_clusters*2); %(x1,y1,x2,y2,...)
cluster_sizes_over_time = cell(num_steps, 1);
attached_phages = false(num_phages, 1);
n_phages_vs_time = zeros(num_steps);

%num_phages = length(phages);
num_bacteria_total = 50;
max_num_clusters = 15;  
max_bacteria_per_cluster = 10; 
phage_positions_over_time = zeros(num_steps, 2 * num_phages);
bact_positions_over_time = zeros(num_steps, 2 * num_bacteria_total);
cluster_positions_over_time = NaN(num_steps, 2 * max_num_clusters);
bact_in_cluster_positions_over_time = NaN(num_steps, 2*num_bacteria_total);
%% Time-stepping loop
disp('----> Begin time iterations')
tic
for k = 1:num_steps

    disp('Time = '); disp(k*dt);

    disp('-------> Compute phages-bacteria interaction forces')
    [bacteriumInteractionForces, phageInteractionForces] = compute_bacteriaPhageInteractionForces(bacteria, phages, epsilon);

    disp('-------> Compute attachments')
    last_time_step = (k == num_steps); 
    [phages, bacteria, attached_phages, n_phages_attached] = compute_attachments(phages, bacteria, clusters, ...
        max_num_bacteria, attached_phages, d_enc1, d_enc2, d_enc3, last_time_step, outputFolder);
    
    disp('-------> Update phage attachments')
    if k == 1
        n_phages_vs_time(k) = n_phages_attached; %attached to bacteria (single or in cluster) only
    else
        n_phages_vs_time(k) = n_phages_attached + n_phages_vs_time(k-1);
    end

    disp('-------> Update phages positions')
    for i = 1:length(phages)
    %parfor i = 1:length(phages)
        if ~attached_phages(i) 
            u_fluid = U_interp([phages(i).position]);
            [phages(i), phageNoiseTerm, phageFluidForce] = phages(i).computeFluidForce(u_fluid);
            phagesTotalForces = phageFluidForce + phageInteractionForces(i);
            phages(i) = phages(i).updateVelocity(phagesTotalForces, dt);
            phages(i) = phages(i).updatePosition(dt, domain, x_min, y_min, x_max, y_max);
        end
    end

    disp('-------> Update bacteria positions')
    for j = 1:length(bacteria) 
    %parfor j = 1:length(bacteria) 
        u_fluid = U_interp([bacteria(j).position]);
        bacteria(j) = bacteria(j).computePropulsionForce();
        bacteria(j) = bacteria(j).computeFluidForce(u_fluid);
        bacteriaTotalForces = bacteria(j).computeTotalBacteriumForce() + bacteriumInteractionForces(j);
        bacteria(j) = bacteria(j).updateVelocity(bacteriaTotalForces, dt);
        bacteria(j) = bacteria(j).updatePosition(dt, domain, x_min, y_min, x_max, y_max);
    end

    disp('-------> Update clusters positions (center of mass)')
    for n = 1:length(clusters)
    %parfor n = 1:length(clusters)
        u_fluid = U_interp([clusters(n).position]);
        clusters(n) = clusters(n).update_friction_coefficient(mu_water, wB);
        clusters(n) = clusters(n).evolveLangevin(u_fluid, dt, kB, T, domain, x_min, y_min, x_max, y_max);
    end
    valid_clusters = arrayfun(@(c) c.size >= 2, clusters);
    cluster_sizes_filtered = arrayfun(@(c) c.size, clusters(valid_clusters))';
    cluster_sizes_over_time{k} = cluster_sizes_filtered(:); %cluster_sizes_over_time{k} = [n₁, n₂, ..., nM], where nᵢ is number of bacteria in cluster i at time step k

    disp('-------> Compute new clusters and bacteria not in cluster')
    bacteria = [bacteria, clusters.bacteria]; % put all bacteria (single and in cluster) in bacteria variable

    [clusters, bacteria] = Cluster.form_clusters(bacteria, threshold, domain, x_min, y_min, x_max, y_max);
    for n = 1:length(clusters)
    %parfor n = 1:length(clusters)
        clusters(n) = clusters(n).update_center_of_mass_and_group_velocity();
        clusters(n).id = n;
    end

    % Print formed clusters updated
    fprintf('At time t = %6f, formed clusters:\n', k*dt);
    for i = 1:length(clusters)
        fprintf('Cluster %d: %d bacteria\n', i, clusters(i).size);
        ids = arrayfun(@(b) b.id, clusters(i).bacteria);
        fprintf('   Bacteria IDs: %s\n', mat2str(ids))
    end

    disp('-------> Plot for real-time visualisation')
    figure(fig1); hold on;

    positionsB_InCluster = [];
    comC = [];
    for i = 1:length(clusters)
        for j = 1:length(clusters(i).bacteria)
            positionsB_InCluster(:, end+1) = clusters(i).bacteria(j).position(:);
        end
        if clusters(i).size > 1
            comC(:, end+1) = clusters(i).position;
        end
    end

    all_bacteria = [clusters.bacteria];
    positionsB_InClu = reshape([all_bacteria.position], 2, []);
    xC = positionsB_InCluster(1, :); %x of bacteria in cluster
    yC = positionsB_InCluster(2, :); %y of bacteria in cluster

    all_clusters = [clusters.position];
    comC = reshape(all_clusters, 2, []);
    xCOM = comC(1, :); %x of com
    yCOM = comC(2, :); %y of com
  
    positionsB = reshape([bacteria.position], 2, []);
    xB = positionsB(1, :); %x of bacteria not in cluster
    yB = positionsB(2, :); %y of bacteria not in cluster

    positionsP = reshape([phages.position], 2, []);
    xP = positionsP(1, :); %x of phages
    yP = positionsP(2, :); %y of phages

    set(h1, 'XData', xC / micron, 'YData', yC / micron );
    set(h2, 'XData', xCOM / micron, 'YData', yCOM / micron);
    set(h3, 'XData', xB / micron, 'YData', yB / micron);
    set(h4, 'XData', xP / micron, 'YData', yP / micron);

    current_time = (k-1) * dt;
    set(time_text1, 'String', sprintf('Time: %.4f s', current_time));
    drawnow; pause(0.2);
    writeVideo(video1, getframe(gcf));

    disp('-------> Save phages and bacteria positions')
%     xP_yP_matrix = [xP; yP];
%     xP_yP_row_fixed_time = reshape(xP_yP_matrix, [2*length(phages),1])';
%     coordP_over_time(k,:) = xP_yP_row_fixed_time;
    
%     xB_yB_matrix = [xB, xC; yB, yC];     %xB_yB_matrix = [xB;yB];
%     xB_yB_row_fixed_time = reshape(xB_yB_matrix, [2*(length(positionsB)+length(all_bacteria)),1])';
%     coordB_over_time(k,:) = xB_yB_row_fixed_time;
    
%     clusters_filtered = clusters([clusters.size] >= 2);
%     num_clusters_k = length(clusters_filtered);
%     xC_yC_matrix = NaN(2, num_clusters_k);  % [2 x num_clusters_k]
%     for c = 1:num_clusters_k
%         xC_yC_matrix(:, c) = clusters_filtered(c).position(:);
%     end
%     xC_yC_row_fixed_time = reshape(xC_yC_matrix, [2*num_clusters_k, 1])';
%     coordC_over_time(k, 1:length(xC_yC_row_fixed_time)) = xC_yC_row_fixed_time;

    phage_positions_over_time = save_phage_positions(phages, k, phage_positions_over_time);
    bact_positions_over_time = save_bacteria_positions(bacteria, k, bact_positions_over_time);
    cluster_positions_over_time = save_cluster_positions(clusters, k, cluster_positions_over_time);
    bact_in_cluster_positions_over_time = save_bacteriaInCluster_positions(clusters, k, ...
        num_bacteria_total, bact_in_cluster_positions_over_time);
end
elapsed_time = toc; % end timer and get elapsed time
fprintf('Total simulation time: %.2f seconds\n', elapsed_time);
%% Post-processing
disp('----> Save data and plot results')
time_vec = (0:num_steps-1) * dt;
%plot_StokesBrinkmanFlow(U_interp, x_min, x_max, y_min, y_max, micron, Omega_X, Omega_Y);
%plot_trajectories_2D(coordP_over_time, coordB_over_time, coordC_over_time, num_phages, num_bacteria, num_clusters_k, micron, Omega_X, Omega_Y);
plot_trajectories_2D(phage_positions_over_time, bact_positions_over_time, cluster_positions_over_time, ...
    bact_in_cluster_positions_over_time, num_phages, num_bacteria, max_num_clusters, micron, Omega_X, Omega_Y);
%save_trajectories(coordP_over_time, coordB_over_time, coordC_over_time, dt, num_steps, num_phages, num_bacteria, num_clusters_k, micron, outputFolder);
save_phage_attachments(n_phages_vs_time, num_steps, dt, outputFolder);

save(fullfile(outputFolder, 'phage_positions_over_time.txt'), 'phage_positions_over_time', '-ascii');
save(fullfile(outputFolder, 'bact_positions_over_time.txt'), 'bact_positions_over_time', '-ascii');
save(fullfile(outputFolder, 'cluster_positions_over_time.txt'), 'cluster_positions_over_time', '-ascii');
save(fullfile(outputFolder, 'bact_in_cluster_positions_over_time.txt'), 'bact_in_cluster_positions_over_time', '-ascii');
%save(fullfile(outputFolder, 'time_vec.txt'), 'time_vec', '-ascii');

%plot_all_trajectories(phage_positions_over_time, bact_positions_over_time, cluster_positions_over_time, time_vec, outputFolder);

close(video1);
disp('DONE')
