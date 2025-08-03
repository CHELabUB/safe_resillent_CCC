clear; close all; clc

%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%
% search detail and running mesh
dt = 0.1; 
beta_range = 0:0.1:2; 
% running details
num_lead_vehicles = 3; 
alpha = 0.4; 
kappa = 0.6; 
sigma = 0.0; 
% modules folder
modules_folder = fullfile(pwd, "modules");
addpath(modules_folder);

profile_names = {"profile1", "profile2", "profile3", "profile4", "profile5", "profile6"};
N = length(profile_names);
source_data = fullfile(fileparts(pwd), "source_data");
for i = 1:N
    Data(i) = load(fullfile(source_data, "8_vehicles_" + profile_names{i} + ".mat")); % load each data file
end

result_root = "perturbation_eqpt_results";
if ~exist(result_root, 'dir')
    mkdir(result_root);
end

run_name = "spectrum_" + num2str(num_lead_vehicles) + "_leading_vehicles";
result_dir = fullfile(result_root, run_name);

run_info.data = profile_names;
run_info.dt = dt;
run_info.beta_range = beta_range;
run_info.num_lead_vehicles = num_lead_vehicles;

n = length(beta_range);
spectrum_opt_data_summary = cell(N,1);
[ii, jj, kk] = ndgrid(1:n, 1:n, 1:n);
beta_combinations = [beta_range(ii(:)); beta_range(jj(:)); beta_range(kk(:))]; % 3 by n^3
beta_combinations = beta_combinations'; % n^3 by 3
nn = size(beta_combinations, 1); % n^3

for i = 1:N
    data = Data(i);
    result_conclusion = zeros(nn, 4);
    run_folder = fullfile(result_dir, profile_names{i});
    tic
    disp('running')
    if ~exist(run_folder, 'dir')
        mkdir(run_folder);
    end
    parfor ind = 1:nn
        beta = beta_combinations(ind, :);
        [J, profile, param] = spectrum_simulation(data, dt, num_lead_vehicles, alpha, beta, kappa, sigma);
        result_conclusion(ind, :) = [beta , J]; % beta 1, beta 2, beta 3, J
        run = struct();
        run.beta = beta;
        run.J = J;
        run.profile = profile;
        run.param = param;
        target_path = fullfile(run_folder, sprintf('%04d.mat', ind));
        save(target_path,"-fromstruct", run);
    end
    toc
    data_summary = fullfile(run_folder, "data_summary.mat");
    save(data_summary, 'result_conclusion');
    spectrum_opt_data_summary{i} = result_conclusion;
end

run_info.spectrum_opt_data_summary = spectrum_opt_data_summary;
summary_path = fullfile(result_dir, 'run_info.mat');
save(summary_path, 'run_info');




