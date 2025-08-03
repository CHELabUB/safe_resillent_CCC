clear; close all; clc

% result directory, all results should be within here

result_root = fullfile(pwd, 'results');
if ~exist(result_root, 'dir')
    mkdir(result_root);
end

% please note that the source data files will be available upon reasonable request, 
% please contact the author for more information.
% Here you can use your on data folder to replace the source data directory to use our simulation pipeline. 

profile_names = {"profile1", "profile2", "profile3", "profile4", "profile5", "profile6"}; 
N = length(profile_names);
source_data = fullfile(fileparts(pwd), "source_data");
modules = addpath(fullfile(pwd,"modules"));
for i = 1:N
    Data(i) = load(fullfile(source_data, "8_vehicles_" + profile_names{i} + ".mat")); % load each data file
end

dt = 0.1; 
beta_range = 0:0.1:2;
use_filter = 0;
vehicle_type = "car";
num_of_lead_vehicles = 3; % number of lead vehicles in the data set
a_under = 4; % underbar of the ego vehicle deceleration capability, change to check the filter effect

run_name = "CCC_" + vehicle_type + "_filter_" + num2str(use_filter) + "_" + num2str(num_of_lead_vehicles) + "_lead_vehicles" + "_a_under_" + num2str(a_under); % name of this run
result_dir = fullfile(result_root, run_name);

% details of this run
run_info.data = profile_names;
run_info.dt = dt;
run_info.beta_range = beta_range;
run_info.use_filter = use_filter;
run_info.vehicle_type = vehicle_type;
run_info.num_of_lead_vehicles = num_of_lead_vehicles;

n = length(beta_range);
CCC_opt_data_summary = cell(N,1);
[ii, jj, kk] = ndgrid(1:n, 1:n, 1:n);
beta_combinations = [beta_range(ii(:)); beta_range(jj(:)); beta_range(kk(:))]; % 3 by n^3
beta_combinations = beta_combinations'; % n^3 by 3
nn = size(beta_combinations, 1); % n^3

% set up workers for further parallel computing
if isempty(gcp('nocreate'))
    numWorkers = min(48, feature('numcores')); % Use up to 24 workers or the number of physical cores
    parpool('local', numWorkers);
else
    delete(gcp('nocreate'));
    numWorkers = min(48, feature('numcores')); % Use up to 24 workers or the number of physical cores
    parpool('local', numWorkers);
end
for i = 1:N
    data = Data(i);
    % CCC_opt_data = cell(n, n, n);
    result_collection = zeros(nn, 6); % [beta1, beta2, beta3, energy_total, unsafe, collision]
    run_folder = fullfile(result_dir, profile_names{i}); % create a folder for each data set
    tic
    if ~exist(run_folder, 'dir')
        mkdir(run_folder)
    end

    parfor ind = 1:nn
        beta = beta_combinations(ind, :); % 1 by 3
        [energy_total, unsafe, collision, profile, param] = v_sim(dt, beta, data, use_filter, vehicle_type, num_of_lead_vehicles, a_under);
        result_collection(ind,:) = [beta, energy_total, unsafe, collision]; % store the result
        run = struct();
        run.beta = beta;
        run.energy = energy_total;
        run.unsafe = unsafe;
        run.collision = collision;
        run.profile = profile;
        run.param = param;
        target_path = fullfile(run_folder, sprintf('%04d.mat', ind)); 
        save(target_path,"-fromstruct", run);
    end
    toc
    % result saving for each run 
    data_summary = fullfile(run_folder, "data_summary.mat");
    save(data_summary, 'result_collection') 
    CCC_opt_data_summary{i} = result_collection; 
end

run_info.CCC_opt_data_summary = CCC_opt_data_summary; 
summary_file = fullfile(result_dir, "run_info.mat");
save(summary_file, 'run_info') % save the summary of this run
