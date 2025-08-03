clear; close all; clc 
% this file is written to extract all the optimized beta combinations of all data folders
% Please note those data folders are not included in the repository, you need to run the simulation pipeline to generate them.
% You need to run the brute force search code in energy_simulation folder twice 
% to get the folder CCC_car_filter_0_3_lead_vehicles_a_under_4 and CCC_car_filter_1_3_lead_vehicles_a_under_4.
% Then run the spectrum analysis code in spectrum_analysis folder.
% the results will be saved in the spectrum_analysis/perturbation_eqpt_results/spectrum_3_leading_vehicles folder.
% If you have any questions, please contact the author. Original data files are available upon reasonable request.


spectrum_run = load(fullfile(fileparts(pwd),'spectrum_analysis','perturbation_eqpt_results','spectrum_3_leading_vehicles', 'run_info.mat'));
spectrum_run_info = spectrum_run.run_info;
spectrum_data_summary = spectrum_run.run_info.spectrum_opt_data_summary;
nn = size(spectrum_run_info.beta_range,2);
N = size(spectrum_data_summary,1);

spectrum_ACC_beta = zeros(N,3);
spectrum_CCC_beta = zeros(N,3);
for i = 1:N
    [min_J, min_idx] = min(spectrum_data_summary{i}(1:nn,4));
    spectrum_ACC_beta(i,:) = spectrum_data_summary{i}(min_idx,1:3);
end
for i = 1:N
    [min_J, min_idx] = min(spectrum_data_summary{i}(:,4));
    spectrum_CCC_beta(i,:) = spectrum_data_summary{i}(min_idx,1:3);
end

mkdir('data');
if ~exist('data', 'dir')
    mkdir('data');
end

save(fullfile('data','spectrum_ACC_beta.mat'), "spectrum_ACC_beta")
save(fullfile('data','spectrum_CCC_beta.mat'), "spectrum_CCC_beta")

filter_0_run_info =  load(fullfile(fileparts(pwd), "energy_simulation", "results", "CCC_car_filter_0_3_lead_vehicles_a_under_4", "run_info.mat"));
filter_1_run_info =  load(fullfile(fileparts(pwd), "energy_simulation", "results", "CCC_car_filter_1_3_lead_vehicles_a_under_4", "run_info.mat"));

filter_1_data_summary = filter_1_run_info.run_info.CCC_opt_data_summary;
filter_0_data_summary = filter_0_run_info.run_info.CCC_opt_data_summary;

ACC_idx = zeros(N,1);
CCC_idx = zeros(N,1);
ACC_data_dir = fullfile('data', 'ACC_data');
if ~exist(ACC_data_dir, 'dir')
    mkdir(ACC_data_dir);
end

CCC_data_dir = fullfile('data', 'CCC_data');
if ~exist(CCC_data_dir, 'dir')
    mkdir(CCC_data_dir);
end
ACC_1 = fullfile(ACC_data_dir, 'filter_1');
if ~exist(ACC_1, 'dir')
    mkdir(ACC_1);
end
ACC_0 = fullfile(ACC_data_dir, 'filter_0');
if ~exist(ACC_0, 'dir')
    mkdir(ACC_0);
end
CCC_1 = fullfile(CCC_data_dir, 'filter_1');
if ~exist(CCC_1, 'dir')
    mkdir(CCC_1);
end
CCC_0 = fullfile(CCC_data_dir, 'filter_0');
if ~exist(CCC_0, 'dir')
    mkdir(CCC_0);
end
CCC_data_dir = fullfile('data', 'CCC_data');
if ~exist(CCC_data_dir, 'dir')
    mkdir(CCC_data_dir);
end
for i = 1:N
    summary = filter_1_data_summary{i};
    ACC_idx(i) = find(summary(:,1) == spectrum_ACC_beta(i,1) & summary(:,2)... 
    == spectrum_ACC_beta(i,2) & summary(:,3) == spectrum_ACC_beta(i,3));
    CCC_idx(i) = find(summary(:,1) == spectrum_CCC_beta(i,1) & summary(:,2)...
    == spectrum_CCC_beta(i,2) & summary(:,3) == spectrum_CCC_beta(i,3));
end

time_mark = {"114049", "115941", "121029", "122850", "124525", "130339"};
for i = 1:N 
    ACC_beta_idx = ACC_idx(i);
    CCC_beta_idx = CCC_idx(i);
    ACC_data_0 = fullfile(fileparts(pwd), "energy_simulation", "results", 'CCC_car_filter_0_3_lead_vehicles_a_under_4', time_mark{i}, sprintf('%04d.mat', ACC_beta_idx));
    ACC_data_1 = fullfile(fileparts(pwd), "energy_simulation", "results", 'CCC_car_filter_1_3_lead_vehicles_a_under_4', time_mark{i}, sprintf('%04d.mat', ACC_beta_idx));
    CCC_data_0 = fullfile(fileparts(pwd), "energy_simulation", "results", 'CCC_car_filter_0_3_lead_vehicles_a_under_4', time_mark{i}, sprintf('%04d.mat', CCC_beta_idx));
    CCC_data_1 = fullfile(fileparts(pwd), "energy_simulation", "results", 'CCC_car_filter_1_3_lead_vehicles_a_under_4', time_mark{i}, sprintf('%04d.mat', CCC_beta_idx));
    ACC_data_0 = load(ACC_data_0);
    ACC_data_1 = load(ACC_data_1);
    CCC_data_0 = load(CCC_data_0);
    CCC_data_1 = load(CCC_data_1);
    ACC_0_save = fullfile(ACC_0, sprintf('%01d.mat', i));
    ACC_1_save = fullfile(ACC_1, sprintf('%01d.mat', i));
    CCC_0_save = fullfile(CCC_0, sprintf('%01d.mat', i));
    CCC_1_save = fullfile(CCC_1, sprintf('%01d.mat', i));
    save(ACC_0_save, 'ACC_data_0');
    save(ACC_1_save, 'ACC_data_1');     
    save(CCC_0_save, 'CCC_data_0');
    save(CCC_1_save, 'CCC_data_1');
end
