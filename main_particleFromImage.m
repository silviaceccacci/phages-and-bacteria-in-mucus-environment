clear all; close all; clc;
%% Add paths
addpath('./particles')
addpath('./fluid')
addpath('./fluid/solverNS/')
addpath('./interpolation')
addpath('./output')
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

load([case_name '_uInterp_1e-4'])

Omega_X = x_max - x_min;
Omega_Y = y_max - y_min;
domain = [Omega_X, Omega_Y];
%% Numerical parameters
disp('----> Read numerical and physical parameters')
dt = 1e-4;
num_steps = 20; %num_steps=60000 to obtain results for t=6sec
%% Physical parameters
mu_water = 10^(-3);     % Dynamic viscosity (Pa s)
rho_water = 10^3;       % Mass density of phage (kg/m^3)
kB = 1.38e-23;          % Boltzmann constant (J/K)
T = 300;                % Temperature (Kelvin)
%% Phages
disp('----> Intialise phages')
rP = 100 * nano;                % Radius of phage (m)
rhoP = 10^5*rho_water;          % Mass density of phage (kg/m^3)
d_enc1 = rP;                    % Encounter distance for pha-bac attachemnt (lyse)
d_enc2 = rP/2;                  % Encounter distance for pha-bacInCl attachemnt (lyse)
d_enc3 = rP/2;                  % Encounter distance for pha-COMcl attachemnt (lyse)
num_phages = 400;

phages_pos = generate_periodic_grid(20, 20, Omega_X, Omega_Y);
phages_pos = phages_pos(1:400,:);

for i = 1:num_phages
    phages(i) = Phage(rP, rhoP, mu_water, kB, T, domain, x_min, y_min, x_max, y_max, phages_pos(i,:));
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

bact_pos = generate_periodic_grid(5, 3, Omega_X, Omega_Y);  
bact_pos = bact_pos(1:num_bacteria,:);

for i = 1:num_bacteria
    bacteria(i) = Bacterium(lB, wB, rhoB, mu_water, kB, T, dt, vB, omega_T, domain, ...
        x_min, y_min, x_max, y_max, bact_pos(i,:));
    bacteria(i).id = i;
    bacteria(i).phages_ids = [];
end
%% Clusters
disp('----> Initialise clusters at initial configuration')
D = 1e-12;
threshold = compute_clustering_threshold(U_interp, x_min, x_max, y_min, y_max, lB, wB, dt, D);
[clusters, bacteria] = Cluster.form_clusters(bacteria, threshold, domain, x_min, y_min, x_max, y_max);

% Print formed clusters
fprintf('Formed %d clusters\n', length(clusters));
for i = 1:length(clusters)
    fprintf('Cluster %d: %d bacteria\n', i, clusters(i).size);
    ids = arrayfun(@(b) b.id, clusters(i).bacteria);
    fprintf('   Bacteria IDs: %s\n', mat2str(ids))
end

for i = 1:length(clusters)
    clusters(i) = clusters(i).update_center_of_mass_and_group_velocity(domain);
    clusters(i).id = i;
end
%% Allocation of variables
disp('----> Allocate variables')
coordP_over_time = zeros(num_steps, num_phages*2);
coordB_over_time = zeros(num_steps, max_num_bacteria*2);
max_clusters = length(bacteria);
coordC_over_time = zeros(num_steps, max_clusters*2); 
cluster_sizes_over_time = cell(num_steps, 1);
attached_phages = false(num_phages, 1);
n_phages_vs_time = zeros(num_steps, 1);

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
phages_radius = cat(1, phages(:).radius);
for k = 1:num_steps

    disp('Time = '); disp(k*dt);

    disp('-------> Compute phages-bacteria interaction forces')
    [bacteriumInteractionForces, phageInteractionForces] = compute_bacteriaPhageInteractionForces(bacteria, phages, epsilon, phages_radius);

    disp('-------> Compute attachments')
    last_time_step = (k == num_steps);
    [phages, bacteria, attached_phages, n_phages_attached] = compute_attachments(phages, bacteria, clusters, ...
        max_num_bacteria, attached_phages, d_enc1, d_enc2, d_enc3, last_time_step, outputFolder);

    disp('-------> Update phage attachments')
    if k == 1
        n_phages_vs_time(k) = n_phages_attached;
    else
        n_phages_vs_time(k) = n_phages_attached + n_phages_vs_time(k-1);
    end

    disp('-------> Update phages positions')
    for i = 1:length(phages)
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
        u_fluid = U_interp([bacteria(j).position]);
        bacteria(j) = bacteria(j).computePropulsionForce();
        bacteria(j) = bacteria(j).computeFluidForce(u_fluid);
        bacteriaTotalForces = bacteria(j).computeTotalBacteriumForce() + bacteriumInteractionForces(j);
        bacteria(j) = bacteria(j).updateVelocity(bacteriaTotalForces, dt);
        bacteria(j) = bacteria(j).updatePosition(dt, domain, x_min, y_min, x_max, y_max);
    end

    disp('-------> Update clusters positions (center of mass)')
    for n = 1:length(clusters)
        u_fluid = U_interp([clusters(n).position]);
        clusters(n) = clusters(n).update_friction_coefficient(mu_water, wB);
        clusters(n) = clusters(n).evolveLangevin(u_fluid, dt, kB, T, domain, x_min, y_min, x_max, y_max);
    end
    valid_clusters = arrayfun(@(c) c.size >= 2, clusters);
    cluster_sizes_filtered = arrayfun(@(c) c.size, clusters(valid_clusters))';
    cluster_sizes_over_time{k} = cluster_sizes_filtered(:); 

    disp('-------> Compute new clusters and bacteria not in cluster')
    if ~isempty(clusters)
        bacteria = [bacteria, clusters.bacteria]; 
    end

    [clusters, bacteria] = Cluster.form_clusters(bacteria, threshold, domain, x_min, y_min, x_max, y_max);
    for n = 1:length(clusters)
        clusters(n) = clusters(n).update_center_of_mass_and_group_velocity(domain);
        clusters(n).id = n;
    end

    num_clusters_at_time(k) = length(clusters);
    current_sizes = arrayfun(@(c) c.size, clusters);
    cluster_sizes_over_time{k} = current_sizes;

    % Print formed clusters updated
    fprintf('At time t = %6f, formed clusters:\n', k*dt);
    for i = 1:length(clusters)
        fprintf('Cluster %d: %d bacteria\n', i, clusters(i).size);
        ids = arrayfun(@(b) b.id, clusters(i).bacteria);
        fprintf('   Bacteria IDs: %s\n', mat2str(ids))
    end

    disp('-------> Save phages and bacteria positions')
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
save_phage_attachments(n_phages_vs_time, num_steps, dt, outputFolder);
save(fullfile(outputFolder, 'phage_positions_over_time.txt'), 'phage_positions_over_time', '-ascii');
save(fullfile(outputFolder, 'bact_positions_over_time.txt'), 'bact_positions_over_time', '-ascii');
save(fullfile(outputFolder, 'cluster_positions_over_time.txt'), 'cluster_positions_over_time', '-ascii');
save(fullfile(outputFolder, 'bact_in_cluster_positions_over_time.txt'), 'bact_in_cluster_positions_over_time', '-ascii');

disp('DONE')
