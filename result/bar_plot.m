clear; close all; clc 

modules_folder = fullfile(pwd, "modules");
addpath(modules_folder);

%%%%%%%%%%%%%%%%%% extracting the spectrum optimal betas %%%%%%%%%%%%%%%%%%%
% PLease note that to use this plotter, you need to run the brute force search code in energy_simluation folder twice 
% to get the folder CCC_car_filter_0_3_lead_vehicles_a_under_4 and CCC_car_filter_1_3_lead_vehicles_a_under_4.
% Then run the spectrum analysis code in spectrum_analysis folder.
% the results will be saved in the spectrum_analysis/pertubation_eqpt_results/spectrum_3_leading_vehicles folder.
% Please note these data folders are not included in the repository, you need to run the simulation pipeline to generate them.
% with the three data fodler mentioend above, you can run this code to plot the bar plot.
% If you have any questions, please contact the author. Original data files are available upon reasonable request.

root_dir = fileparts(pwd);
spectrum_run = load(fullfile(root_dir, 'spectrum_analysis','pertubation_eqpt_results','spectrum_3_leading_vehicles', 'run_info.mat'));
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

%%%%%%%%%%%%%%%%%% extracting the smoothed ACC and CCC data %%%%%%%%%%%%%%%%%%%
filter_0_run_info =  load(fullfile(root_dir, "energy_simulation", "results", "CCC_car_filter_0_3_lead_vehicles_a_under_4", "run_info.mat"));
filter_1_run_info =  load(fullfile(root_dir, "energy_simulation", "results", "CCC_car_filter_1_3_lead_vehicles_a_under_4", "run_info.mat"));

filter_1_data_summary = filter_1_run_info.run_info.CCC_opt_data_summary;  % only need the ACC data summary with the filter_1
filter_0_data_summary = filter_0_run_info.run_info.CCC_opt_data_summary;

ACC_idx = zeros(N,N);
CCC_idx = zeros(N,N);
ACC_cross_energy = zeros(N,N);
CCC_cross_energy = zeros(N,N);
% ACC_data_dir = fullfile('data', 'ACC_data');

for i = 1:N
    for j = 1:N
    summary = filter_1_data_summary{j};
    ACC_idx(i, j) = find(summary(:,1) == spectrum_ACC_beta(i,1) & summary(:,2)... 
    == spectrum_ACC_beta(i,2) & summary(:,3) == spectrum_ACC_beta(i,3));
    ACC_cross_energy(i, j) = filter_1_data_summary{j}(ACC_idx(i,j),4);
    CCC_idx(i, j) = find(summary(:,1) == spectrum_CCC_beta(i,1) & summary(:,2)...
    == spectrum_CCC_beta(i,2) & summary(:,3) == spectrum_CCC_beta(i,3));
    CCC_cross_energy(i, j) = filter_1_data_summary{j}(CCC_idx(i,j),4);
    end
end

%% extracting brute force search optimals

brute_ACC_beta = zeros(N,3);
brute_CCC_beta = zeros(N,3);
brute_ACC_energy = zeros(N,1);
brute_CCC_energy = zeros(N,1);
brute_ACC_idx = zeros(N,1);
brute_CCC_idx = zeros(N,1);

for i = 1:N
    [min_J, min_idx] = min(filter_1_data_summary{i}(1:nn,4));
    brute_ACC_idx(i) = min_idx;
    brute_ACC_beta(i,:) = filter_1_data_summary{i}(min_idx,1:3);
    brute_ACC_energy(i) = filter_1_data_summary{i}(min_idx,4);
end
for i = 1:N
    [min_J, min_idx] = min(filter_1_data_summary{i}(:,4));
    brute_CCC_idx(i) = min_idx;
    brute_CCC_beta(i,:) = filter_1_data_summary{i}(min_idx,1:3);
    brute_CCC_energy(i) = filter_1_data_summary{i}(min_idx,4);
end

%% bar code
% layer: CCC_brute -> CCC_spectrum -> ACC_brute -> ACC_spectrum
brute_ACC_energy = repmat(brute_ACC_energy', N, 1);
brute_CCC_energy = repmat(brute_CCC_energy', N, 1);
% get the percentage of the improvement
ACC_spectrum_improvement = (ACC_cross_energy - ACC_cross_energy)./ACC_cross_energy .* 100;
ACC_brute_improvement = -(brute_ACC_energy - ACC_cross_energy)./ACC_cross_energy .* 100;
CCC_spectrum_improvement = -(CCC_cross_energy - ACC_cross_energy)./ACC_cross_energy .* 100;
CCC_brute_improvement = -(brute_CCC_energy - ACC_cross_energy)./ACC_cross_energy .* 100;

colors = {
    %[82, 245, 0]./255;  
    [255, 170, 170]./255; 
    [191, 40, 40]./255; 
    [32, 10, 144]./255;
    [204, 229, 255]./255;
    [44, 5, 241]./255
    %[82, 245, 0]./255;  
%rgb(44, 5, 241)
};


CB = brute_CCC_energy./1000;
CS = CCC_cross_energy./1000;
AB  = brute_ACC_energy./1000;
AS = ACC_cross_energy./1000;


fig_energy = figure(1); clf;
set_fig(fig_energy, 16, 3, "energy_improvement_compare")
% figureset(gcf, 'Height',3,'Width', 16)
x = 1:5;
x_tick = x;
bar(x, AS(1, [2:6]),'FaceColor', colors{2}, 'BarWidth', 0.8);
hold on
bar(x + 0.1, AB(1, [2:6]), 'FaceColor', colors{1}, 'BarWidth', 0.6);
bar(x + 0.1, CS(1, [2:6]), 'FaceColor', colors{3}, 'BarWidth', 0.6);
bar(x + 0.2, CB(1, [2:6]), 'FaceColor', colors{4}, 'BarWidth', 0.4);
xticks([1 2 3 4 5]);  % Set specific xticks

x = x + 6;
x_tick = [x_tick, x];
bar(x, AS(2, [1, 3:6]),'FaceColor', colors{2}, 'BarWidth', 0.8);
hold on
bar(x + 0.1, AB(2, [1, 3:6]), 'FaceColor', colors{1}, 'BarWidth', 0.6);
bar(x + 0.1, CS(2, [1, 3:6]), 'FaceColor', colors{3}, 'BarWidth', 0.6);
bar(x + 0.2, CB(2, [1, 3:6]), 'FaceColor', colors{4}, 'BarWidth', 0.4);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, AS(3, [1, 2, 4:6]),'FaceColor', colors{2}, 'BarWidth', 0.8);
hold on
bar(x + 0.1, AB(3, [1, 2, 4:6]), 'FaceColor', colors{1}, 'BarWidth', 0.6);
bar(x + 0.1, CS(3, [1, 2, 4:6]), 'FaceColor', colors{3}, 'BarWidth', 0.6);
bar(x + 0.2, CB(3, [1, 2, 4:6]), 'FaceColor', colors{4}, 'BarWidth', 0.4);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, AS(4, [1:3, 5, 6]),'FaceColor', colors{2}, 'BarWidth', 0.8);
hold on
bar(x + 0.1, AB(4, [1:3, 5, 6]), 'FaceColor', colors{1}, 'BarWidth', 0.6);
bar(x + 0.1, CS(4, [1:3, 5, 6]), 'FaceColor', colors{3}, 'BarWidth', 0.6);
bar(x + 0.2, CB(4, [1:3, 5, 6]), 'FaceColor', colors{4}, 'BarWidth', 0.4);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, AS(5, [1:4, 6]),'FaceColor', colors{2}, 'BarWidth', 0.8);
hold on
bar(x + 0.1, AB(5, [1:4, 6]), 'FaceColor', colors{1}, 'BarWidth', 0.6);
bar(x + 0.1, CS(5, [1:4, 6]), 'FaceColor', colors{3}, 'BarWidth', 0.6);
bar(x + 0.2, CB(5, [1:4, 6]), 'FaceColor', colors{4}, 'BarWidth', 0.4);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, AS(6, [1:5]),'FaceColor', colors{2}, 'BarWidth', 0.8);
hold on
bar(x + 0.1, AB(6, [1:5]), 'FaceColor', colors{1}, 'BarWidth', 0.6);
bar(x + 0.1, CS(6, [1:5]), 'FaceColor', colors{3}, 'BarWidth', 0.6);
bar(x + 0.2, CB(6, [1:5]), 'FaceColor', colors{4}, 'BarWidth', 0.4);
xticks(x_tick);

xticklabels({'2', '3', '4', '5', '6', '1', '3', '4', '5', '6', '1', '2', '4', '5', '6', '1', '2', '3', '5', '6', ...
'1', '2', '3', '4', '6', '1', '2', '3', '4', '5' });  % Set custom labels for xticks

% ylabel('$w [\frac{kJ}{kg}]$')
% xlabel('Set number')
set(gca, 'FontName', 'Times', 'FontSize', 10, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 16);


fig_improve = figure(2); clf
set_fig(fig_improve, 16, 2, "robust_improvement_compare")
% figureset(gcf, 'Height',3,'Width', 16)
percentage = (AS - CS) ./AS .* 100;
x = 1:5;
x_tick = x;
bar(x, percentage(1, [2:6]),'FaceColor', colors{5}, 'BarWidth', 0.8);
hold on
xticks([1 2 3 4 5]);  % Set specific xticks

x = x + 6;
x_tick = [x_tick, x];
bar(x, percentage(2, [1, 3:6]),'FaceColor', colors{5}, 'BarWidth', 0.8);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, percentage(3, [1, 2, 4:6]),'FaceColor', colors{5}, 'BarWidth', 0.8);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, percentage(4, [1:3, 5, 6]),'FaceColor', colors{5}, 'BarWidth', 0.8);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, percentage(5, [1:4, 6]),'FaceColor', colors{5}, 'BarWidth', 0.8);
xticks(x_tick);

x = x + 6;
x_tick = [x_tick, x];
bar(x, percentage(6, [1:5]),'FaceColor', colors{5}, 'BarWidth', 0.8);
xticks(x_tick);

xticklabels({'2', '3', '4', '5', '6', '1', '3', '4', '5', '6', '1', '2', '4', '5', '6', '1', '2', '3', '5', '6', ...
'1', '2', '3', '4', '6', '1', '2', '3', '4', '5' });  % Set custom labels for xticks

% ylabel('Energy Improvement$\%$')
% xlabel('Set number')
set(gca, 'FontName', 'Times', 'FontSize', 10, 'TickLabelInterpreter', 'latex');
set(gca, 'FontSize', 16);



full_save_path = @(pic_name) fullfile("paper_fig", pic_name + ".png");
exportgraphics(fig_energy, full_save_path("w_compare"),'BackgroundColor','none','ContentType','vector');
exportgraphics(fig_improve, full_save_path("robust analysis"),'BackgroundColor','none','ContentType','vector');

function plot_bars(CCC_brute, CCC_spectrum, ACC_brute, ACC_spectrum, idx, colors)
    hold on; box on; 
    
    figureset(gcf, 'Height',4,'Width',4);
    switch idx
        case 1
            s_idx = 2:6;
        case 2
            s_idx = [1, 3:6];
        case 3
            s_idx = [1, 2, 4:6];
        case 4
            s_idx = [1:3, 5, 6];
        case 5
            s_idx = [1:4, 6];
        case 6
            s_idx = 1:5;
    end
    CB = CCC_brute(idx, s_idx)./1000;
    CS = CCC_spectrum(idx, s_idx)./1000;
    AB = ACC_brute(idx, s_idx)./1000;
    AS = ACC_spectrum(idx, s_idx)./1000;
    
    bar_to_plot = [CB;CS-CB; AB-CS;AS-AB];
    h = bar(bar_to_plot', 'stacked', 'BarWidth',1);
    
    
    h(1).FaceColor = colors{1};
    h(2).FaceColor = colors{2};
    h(3).FaceColor = colors{3};
    h(4).FaceColor = colors{4};
    % Set the color of each bar
    
    % Add labels, legend, etc. (optional)
    xticks(1:length(s_idx)); % Set the x-ticks to match the number of bars
    xticklabels(arrayfun(@num2str, s_idx, 'UniformOutput', false)); % Set x-tick labels to s_idx values


    % Add labels, legend, etc.
    xlabel('Set number', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel('$\omega$', 'Interpreter', 'latex', 'FontSize', 12);
    title(['Energy comparison for set ', num2str(idx)], 'Interpreter', 'latex', 'FontSize', 14);
    legend({'CCC Brute', 'CCC Spectrum', 'ACC Brute', 'ACC Spectrum'}, ...
        'Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);

    % Set global font properties for the axes
    set(gca, 'FontName', 'Times', 'FontSize', 10, 'TickLabelInterpreter', 'latex');
    set(gca, 'FontSize', 20);
end
