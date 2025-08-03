close; clear; clc 
modules_folder = fullfile(pwd, "modules");
addpath(modules_folder);

set(0,'defaultLineLinewidth', 2); % this way you don't need to do linewith = num_col everywhere
set(0,'defaultAxesFontSize', 12); % this way you don't need to do fontsize everywhere
set(groot, 'DefaultTextInterpreter', 'latex');
set(groot, 'DefaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%%
full_save_path = @(pic_name) fullfile("paper_fig", pic_name + ".png");
% This script is written to investigate the SF impact
% for the sake of space, this data only contains run details of the ACC and CCC spectrum optimal's run details
% data structure:
% data
%    |spectrum_ACC_beta.mat : recording all optimal ACC betas from spectrum 
%    |spectrum_CCC_beta.mat : recording all optimal CCC betas from spectrum
%    |ACC_data  (ACC's are extracted from the CCC)
%          |filter_0 (all ACC data runs without SF)
%               |1.mat  (optimal spectrum ACC beta's run information from energy brute force without SF)
%               |num_col.mat
%                 :
%               |6.mat
%          |filter_1 (all ACC data runs with SF)
%               |1.mat  (optimal spectrum ACC beta's run information from energy brute force with SF)
%               |num_col.mat
%                 :
%               |6.mat
%    |CCC_data
%          |filter_0 (all CCC data runs without SF)
%               |1.mat  (optimal spectrum CCC beta's run information from energy brute force without SF)
%               |num_col.mat
%                 :
%               |6.mat
%          |filter_1 (all CCC data runs with SF)
%               |1.mat  (optimal spectrum CCC beta's run information from energy brute force with SF)
%               |num_col.mat
%                 :
%               |6.mat

% user input
run_interest = 1; 
% which dataset interested to compare, decide to use 1 because it is the one that CCC got hurt with filter while ACC got help.

% loading data, within data, it contains all running details within data (including params and profile etc)
% PLease note that the data folder is not included in the repository, 
% you need to run the file_extract.m to extract the optimal results after you have the brute force search results and spectrum analysis results.

data_root = "data";
fig_name = data_root + "_" + num2str(run_interest);

% ACC data
ACC_filter_0 = fullfile(data_root,'ACC_data','filter_0',strcat(num2str(run_interest),'.mat'));
ACC_filter_1 = fullfile(data_root,'ACC_data','filter_1',strcat(num2str(run_interest),'.mat'));
ACC_filter_0 = load(ACC_filter_0); ACC_filter_0 = ACC_filter_0.ACC_data_0;
ACC_filter_1 = load(ACC_filter_1); ACC_filter_1 = ACC_filter_1.ACC_data_1;
% CCC data
CCC_filter_0 = fullfile(data_root,'CCC_data','filter_0',strcat(num2str(run_interest),'.mat'));
CCC_filter_1 = fullfile(data_root,'CCC_data','filter_1',strcat(num2str(run_interest),'.mat'));      
CCC_filter_0 = load(CCC_filter_0); CCC_filter_0 = CCC_filter_0.CCC_data_0;
CCC_filter_1 = load(CCC_filter_1); CCC_filter_1 = CCC_filter_1.CCC_data_1;

% shift time, do this for visualization purpose
ACC_filter_0.profile.time = ACC_filter_0.profile.time - ACC_filter_0.profile.time(1);
ACC_filter_1.profile.time = ACC_filter_1.profile.time - ACC_filter_1.profile.time(1);
CCC_filter_0.profile.time = CCC_filter_0.profile.time - CCC_filter_0.profile.time(1);
CCC_filter_1.profile.time = CCC_filter_1.profile.time - CCC_filter_1.profile.time(1);
%%
%% Calculate metrics
% profile = CCC_filter_1.profile;
% skeleton, get cumulative profile using building functions
get_cum_profile = @(profile, int_value) cumtrapz(profile.time, int_value(profile));

% total energy, should be the same as the simulation energy profile
w_fun = @(profile) [max(profile.dv + profile.resist, 0) .* profile.V(1:end-1), 0];
get_w = @(run) get_cum_profile(run.profile, w_fun);

w_cum_ccc_with_filter = get_w(CCC_filter_1) / 1000; % J/kg to kJ/kg
w_cum_ccc_no_filter = get_w(CCC_filter_0) / 1000;   % J/kg to kJ/kg
w_cum_acc_with_filter = get_w(ACC_filter_1) / 1000; % J/kg to kJ/kg
w_cum_acc_no_filter = get_w(ACC_filter_0) / 1000;   % J/kg to kJ/kg

% wasted energy, does not counted in the energy calculation, but needed to be recovered
w_brake_energy = @(profile) [max(-(profile.dv + profile.resist), 0) .* profile.V(1:end-1), 0];
get_w_brake = @(run) get_cum_profile(run.profile, w_brake_energy);

w_brake_ccc_with_filter = get_w_brake(CCC_filter_1) / 1000; % J/kg to kJ/kg
w_brake_ccc_no_filter = get_w_brake(CCC_filter_0) / 1000; % J/kg to kJ/kg
w_brake_acc_with_filter = get_w_brake(ACC_filter_1) / 1000; % J/kg to kJ/kg
w_brake_acc_no_filter = get_w_brake(ACC_filter_0) / 1000; % J/kg to kJ/kg

% energy profile when acc kicked in, when a_ccc > a_cbf, the SF kicks in
w_filter_energy = @(profile) [max(profile.dv + profile.resist, 0) .* profile.V(1:end-1) .* ...
                              max(profile.a_ccc - profile.a_cbf, 0), 0];
get_w_filter_energy = @(run) get_cum_profile(run.profile, w_filter_energy);

w_filter_energy_ccc_with_filter = get_w_filter_energy(CCC_filter_1) / 1000; % J/kg to kJ/kg
w_filter_energy_acc_with_filter = get_w_filter_energy(ACC_filter_1) / 1000; % J/kg to kJ/kg

% average kick in energy
% a_ccc - profile.a_safe should be non negative
a_mean_kick_in = @(profile) [profile.a_ccc - profile.a_safe, 0];
a_mean_kick_in_positive = @(profile) [max((profile.a_ccc - profile.a_safe) .* (profile.dv + profile.resist>=0), 0), 0];
a_mean_kick_in_negative = @(profile) [max((profile.a_ccc - profile.a_safe) .* (profile.dv + profile.resist<0), 0), 0];

get_a_mean_kick_in = @(run) get_cum_profile(run.profile, a_mean_kick_in) / (run.profile.time(end) - run.profile.time(1)); % m/s^2
get_a_mean_kick_in_positive = @(run) get_cum_profile(run.profile, a_mean_kick_in_positive) / (run.profile.time(end) - run.profile.time(1));
get_a_mean_kick_in_negative = @(run) get_cum_profile(run.profile, a_mean_kick_in_negative) / (run.profile.time(end) - run.profile.time(1));

a_mean_kick_in_ccc = get_a_mean_kick_in(CCC_filter_1); % m/s^2
a_mean_kick_in_ccc_positive = get_a_mean_kick_in_positive(CCC_filter_1); % m/s^2
a_mean_kick_in_ccc_negative = get_a_mean_kick_in_negative(CCC_filter_1); % m/s^2

a_mean_kick_in_acc = get_a_mean_kick_in(ACC_filter_1); % m/s^2
a_mean_kick_in_acc_positive = get_a_mean_kick_in_positive(ACC_filter_1); % m/s^2
a_mean_kick_in_acc_negative = get_a_mean_kick_in_negative(ACC_filter_1); % m/s^2

%% Unsafe percentage
% need a tolerance to determine if the profile is safe or not, this is because of the time delay
% threshold is necessary to claim the guarantee, the gap is due to the delay.
threshold = 0.65; % need 0.65 for run 1 and 2, 0.5 for the rest is good enough

max(ACC_filter_1.profile.b_hat - ACC_filter_1.profile.H(1:end-1))
max(CCC_filter_1.profile.b_hat - CCC_filter_1.profile.H(1:end-1))

% threshold = 0.0;
h_not_safe = @(profile) [profile.b_hat - profile.H(1:end-1) >= threshold, 0]; % h_not_safe is the distance to the safe distance
h_not_safe_margin = @(profile) [max(profile.b_hat - profile.H(1:end-1) - threshold, 0), 0];
% h_not_safe = @(profile) [profile.b_hat > profile.H, 0]; % h_not_safe is the distance to the safe distance
get_h_not_safe = @(run) get_cum_profile(run.profile, h_not_safe) / (run.profile.time(end) - run.profile.time(1)); % [%]
get_h_not_safe_margin = @(run) get_cum_profile(run.profile, h_not_safe_margin) / (run.profile.time(end) - run.profile.time(1)); % [m]

h_not_safe_ccc_with_filter = get_h_not_safe(CCC_filter_1); 
h_not_safe_ccc_no_filter = get_h_not_safe(CCC_filter_0);  
h_not_safe_acc_with_filter = get_h_not_safe(ACC_filter_1); 
h_not_safe_acc_no_filter = get_h_not_safe(ACC_filter_0); 

h_not_safe_ccc_with_filter_margin = get_h_not_safe_margin(CCC_filter_1); 
h_not_safe_ccc_no_filter_margin = get_h_not_safe_margin(CCC_filter_0); 
h_not_safe_acc_with_filter_margin = get_h_not_safe_margin(ACC_filter_1); 
h_not_safe_acc_no_filter_margin = get_h_not_safe_margin(ACC_filter_0); 

%%
% to make plotter, we need to extract the profiles
time = ACC_filter_0.profile.time; % all profiles share the same time

% four panels, v, D, u, w, extend more if needed, num_row by num_col left is ACC, right is CCC

fig = figure(1);
set_fig(fig, 18, 7, num2str(run_interest) + "ACC vs CCC, unsafe vs safe");
% set(fig, "Name", num2str(run_interest) + "ACC vs CCC, unsafe vs safe")
% set(gcf, "unit", "inches");
% ps = get(gcf, "Position");
% width = 18;
% height = 7;
% set(gcf, "Position", [ps(1)/7, ps(2)/7, width, height])

num_col = 2; % number of columns
num_row = 3; % number of rows


% allow quick changes of color and line style
acc_no_filter_color = 'k'; % ACC no SF color
acc_filter_color = 'r:'; % ACC SF color
acc_filter_pre_sf_color = 'k'; % ACC SF color before SF kicks in
acc_filter_color_2 = 'k'; % ACC SF color additional

ccc_no_filter_color = 'k'; % CCC no SF color
ccc_filter_color = 'b-.'; % CCC SF color
ccc_filter_pre_sf_color = 'k'; % ACC SF color before SF kicks in
ccc_filter_color_2 = 'k'; % CCC SF color additional

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ACC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(num_row,num_col,1); hold on; grid on; box on;
plot(time,ACC_filter_0.profile.V, acc_no_filter_color);
plot(time,ACC_filter_1.profile.V, acc_filter_color);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('ACC w/o SF','ACC w/ SF');

subplot(num_row,num_col,3); hold on; grid on; box on;
plot(time,ACC_filter_0.profile.H, acc_no_filter_color);
plot(time,ACC_filter_1.profile.H, acc_filter_color);
% plot(time(1:end-1),ACC_filter_0.profile.b_hat, acc_filter_color_2, "LineWidth", 1.5)
xlabel('Time (s)'); ylabel('D (m)');
legend('ACC w/o SF','ACC w/ SF');
% legend('ACC w/o SF','ACC w/ SF', "B for ACC w/o SF");

subplot(num_row,num_col,5); hold on; grid on; box on;
plot(time(1:end-1), ACC_filter_1.profile.a_ccc, acc_filter_pre_sf_color);
plot(time(1:end-1), ACC_filter_1.profile.a_safe, acc_filter_color);


xlabel('Time (s)'); ylabel('a (m/s^2)');
legend('ACC (before SF)', 'ACC a safe','Location','northwest');

% subplot(num_row,num_col,7); hold on; grid on; box on;
% plot(time, h_not_safe_acc_no_filter, acc_no_filter_color);
% plot(time, h_not_safe_acc_with_filter, acc_filter_color);
% xlabel('Time (s)'); ylabel('unsafe');
% legend('ACC w/o SF','ACC w/ SF', 'Location','northwest');

% subplot(num_row,num_col,9); hold on; grid on; box on;
% plot(time, h_not_safe_acc_no_filter_margin, acc_no_filter_color);
% plot(time, h_not_safe_acc_with_filter_margin, acc_filter_color);
% xlabel('Time (s)'); ylabel('unsafe margin [m]');
% legend('ACC w/o SF','ACC w/ SF', 'Location','northwest');

% plot(time,ACC_filter_0.profile.energy/1000, acc_no_filter_color);
% plot(time,ACC_filter_1.profile.energy/1000, acc_filter_color);
% % plot(time, w_cum_acc_no_filter,"c--");
% % plot(time, w_cum_acc_with_filter,"m--");
% xlabel('Time (s)'); ylabel('w (kJ/kg)');
% legend('ACC w/o SF','ACC w/ SF', 'Location','northwest');
% 
% subplot(num_row,num_col,9); hold on; grid on; box on;
% plot(time, w_brake_acc_no_filter, acc_no_filter_color);
% plot(time, w_brake_acc_with_filter, acc_filter_color);
% xlabel('Time (s)'); ylabel('w brake (kJ/kg)');
% legend('ACC w/o SF','ACC w/ SF', 'Location','northwest');

% subplot(num_row,num_col,11); hold on; grid on; box on;
% plot(time, w_filter_energy_acc_with_filter, acc_filter_color);
% xlabel('Time (s)'); ylabel('w SF (kJ/kg)');
% legend('ACC w/ SF', 'Location','northwest');

% subplot(num_row,num_col,13); hold on; grid on; box on;
% plot(time(1:end-1),ACC_filter_0.profile.u, acc_no_filter_color);
% plot(time(1:end-1),ACC_filter_1.profile.u, acc_filter_color);
% xlabel('Time (s)'); ylabel('u (m/s)');
% legend('ACC w/o SF','ACC w/ SF');


% subplot(num_row,num_col,15); hold on; grid on; box on;
% plot(time, a_mean_kick_in_acc, acc_filter_color);
% plot(time, a_mean_kick_in_acc_positive, "g--");
% plot(time, a_mean_kick_in_acc_negative, "k--");
% xlabel('Time (s)'); ylabel('a (m/s^2)');
% legend('ACC a mean kick in','ACC a mean kick in positive','ACC a mean kick in negative', 'Location','northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CCC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plot_idx = 2;
subplot(num_row,num_col,2); hold on; grid on; box on;
plot(time,CCC_filter_0.profile.V,ccc_no_filter_color);
plot(time,CCC_filter_1.profile.V,ccc_filter_color);
xlabel('Time (s)'); ylabel('Velocity (m/s)');
legend('CCC w/o SF','CCC w/ SF');

subplot(num_row,num_col,4); hold on; grid on; box on;
plot(time,CCC_filter_0.profile.H,ccc_no_filter_color);                        
plot(time,CCC_filter_1.profile.H,ccc_filter_color);
% plot(time(1:end-1),CCC_filter_0.profile.b_hat, ccc_filter_color_2, 'LineWidth', 1.5)
xlabel('Time (s)'); ylabel('D (m)');
legend('CCC w/o SF','CCC w/ SF');
% legend('CCC w/o SF','CCC w/ SF', "B for CCC w/o SF");
subplot(num_row,num_col,6); hold on; grid on; box on;
plot(time(1:end-1), CCC_filter_1.profile.a_ccc, ccc_filter_pre_sf_color);
plot(time(1:end-1), CCC_filter_1.profile.a_safe, ccc_filter_color);
xlabel('Time (s)'); ylabel('a_d (m/s^2)');
legend('CCC (before SF)', 'CCC a safe','Location','northwest');

% subplot(num_row,num_col,8); hold on; grid on; box on;
% plot(time, h_not_safe_ccc_no_filter, ccc_no_filter_color);
% plot(time, h_not_safe_ccc_with_filter, ccc_filter_color);
% xlabel('Time (s)'); ylabel('unsafe');
% legend('CCC w/o SF','CCC w/ SF', 'Location','northwest');

% subplot(num_row,num_col,10); hold on; grid on; box on;
% plot(time, h_not_safe_ccc_no_filter_margin, ccc_no_filter_color);
% plot(time, h_not_safe_ccc_with_filter_margin, ccc_filter_color);
% xlabel('Time (s)'); ylabel('unsafe margin [m]');
% legend('CCC w/o SF','CCC w/ SF', 'Location','northwest');

% plot(time,CCC_filter_0.profile.energy/1000,ccc_no_filter_color);
% plot(time,CCC_filter_1.profile.energy/1000,ccc_filter_color);
% % plot(time, w_cum_ccc_no_filter,"c--");
% % plot(time, w_cum_ccc_with_filter,"m--");
% xlabel('Time (s)'); ylabel('w (J/kg)');
% legend('CCC w/o SF','CCC w/ SF', 'Location','northwest');
% 
% subplot(num_row,num_col,10); hold on; grid on; box on;
% plot(time, w_brake_ccc_no_filter,ccc_no_filter_color);
% plot(time, w_brake_ccc_with_filter,ccc_filter_color);
% xlabel('Time (s)'); ylabel('w brake (kJ/kg)');
% legend('CCC w/o SF','CCC w/ SF', 'Location','northwest');

% subplot(num_row,num_col,12); hold on; grid on; box on;
% plot(time(1:end-1),CCC_filter_0.profile.u,ccc_no_filter_color);
% plot(time(1:end-1),CCC_filter_1.profile.u,ccc_filter_color);
% xlabel('Time (s)'); ylabel('u (m/s)');
% legend('CCC w/o SF','CCC w/ SF');

% subplot(num_row,num_col,14); hold on; grid on; box on;

% plot(time, w_filter_energy_ccc_with_filter,ccc_filter_color);
% xlabel('Time (s)'); ylabel('w SF (kJ/kg)');
% legend('CCC w/ SF', 'Location','northwest');

% subplot(num_row,num_col,16); hold on; grid on; box on;
% plot(time, a_mean_kick_in_ccc, ccc_filter_color);
% plot(time, a_mean_kick_in_ccc_positive, "g--");
% plot(time, a_mean_kick_in_ccc_negative, "k--");
% xlabel('Time (s)'); ylabel('a (m/s^2)');
% legend('CCC a mean kick in','CCC a mean kick in positive','CCC a mean kick in negative', 'Location','northwest');


%% make consistent y-axis
for idx = 1:num_row
    ylim_all = [inf, -inf];
    for jdx = 1:2
        subplot(num_row,num_col,2* idx -2 + jdx);
        ylim_value = ylim;
        ylim_all(1) = min(ylim_all(1), ylim_value(1));
        ylim_all(2) = max(ylim_all(2), ylim_value(2));
    end
    for jdx = 1:2
        subplot(num_row,num_col,2* idx -2 + jdx);
        ylim(ylim_all);
    end
end

% unify all energy profiles to compare
% energy_idx = (7:10);
% ylim_all = [0, 0];
% for idx = energy_idx
%     subplot(num_row,num_col,idx);
%     % ylim = get(gca,'ylim');
%     ylim_value = ylim;
%     ylim_all(2) = max(ylim_all(2), ylim_value(2));
% end
% ylim_all(2) = ceil(ylim_all(2));
% for idx = energy_idx
%     subplot(num_row,num_col,idx);
%     % ylim(2) = ylim_all(2);
%     ylim(ylim_all);
% end


% Link all x-axes
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

%% for paper fig, plot ACC and CCC panel separately with no legend or xlabel ylabel
ylim_v = [0, 35];
ylim_d = [0, 80];
ylim_a = [-5, 5];

fig_ACC = figure(2);
set_fig(fig_ACC, 8, 7, num2str(run_interest) + "ACC w/o SF vs ACC w/ SF");
subplot(3,1,1); hold on; grid on; box on;
plot(time,ACC_filter_0.profile.H, acc_no_filter_color);
plot(time,ACC_filter_1.profile.H, acc_filter_color);
plot([10, 50], [70, 70], acc_no_filter_color);
plot([250, 290], [70, 70], acc_filter_color);
ylim(ylim_d);
subplot(3,1,2); hold on; grid on; box on;
plot(time,ACC_filter_0.profile.V, acc_no_filter_color);
plot(time,ACC_filter_1.profile.V, acc_filter_color);
plot([10, 50], [5, 5], acc_no_filter_color);
plot([250, 290], [5, 5], acc_filter_color);
ylim(ylim_v);
subplot(3,1,3); hold on; grid on; box on;
plot(time(1:end-1), ACC_filter_1.profile.a_ccc, acc_filter_pre_sf_color);
plot(time(1:end-1), ACC_filter_1.profile.a_safe, acc_filter_color);
plot([10, 50], [3, 3], acc_no_filter_color);
plot([250, 290], [3, 3], acc_filter_color);
ylim(ylim_a);
linkaxes(findall(gcf, 'Type', 'axes'), 'x');

fig_CCC = figure(3);
set_fig(fig_CCC, 8, 7, num2str(run_interest) + "CCC w/o SF vs CCC w/ SF");
subplot(3,1,1); hold on; grid on; box on;
plot(time,CCC_filter_0.profile.H,ccc_no_filter_color);
plot(time,CCC_filter_1.profile.H,ccc_filter_color);
plot([10, 50], [10, 10], ccc_no_filter_color);
plot([250, 290], [10, 10], ccc_filter_color);
ylim(ylim_d);
subplot(3,1,2); hold on; grid on; box on;
plot(time,CCC_filter_0.profile.V,ccc_no_filter_color);
plot(time,CCC_filter_1.profile.V,ccc_filter_color);
plot([10, 50], [5, 5], ccc_no_filter_color);
plot([250, 290], [5, 5], ccc_filter_color);
ylim(ylim_v);
subplot(3,1,3); hold on; grid on; box on;
plot(time(1:end-1), CCC_filter_1.profile.a_ccc, ccc_filter_pre_sf_color);
plot(time(1:end-1), CCC_filter_1.profile.a_safe, ccc_filter_color);
plot([10, 50], [-3.5, -3.5], ccc_no_filter_color);
plot([250, 290], [-3.5, -3.5], ccc_filter_color);
ylim(ylim_a);
linkaxes(findall(gcf, 'Type', 'axes'), 'x');



%%
%% stats
table_row = ["w/o SF"; "w/ SF"];
table_col = ["ACC", "CCC"];
table_h_flag = [h_not_safe_acc_no_filter(end), h_not_safe_ccc_no_filter(end);
               h_not_safe_acc_with_filter(end), h_not_safe_ccc_with_filter(end)];
table_h_margin = [h_not_safe_acc_no_filter_margin(end), h_not_safe_ccc_no_filter_margin(end);
               h_not_safe_acc_with_filter_margin(end), h_not_safe_ccc_with_filter_margin(end)];
table_w = [w_cum_acc_no_filter(end), w_cum_ccc_no_filter(end);
           w_cum_acc_with_filter(end), w_cum_ccc_with_filter(end)];
table_w_brake = [w_brake_acc_no_filter(end), w_brake_ccc_no_filter(end);
                 w_brake_acc_with_filter(end), w_brake_ccc_with_filter(end)];

disp(table_w)
disp(table_w_brake)
disp(table_h_flag)
disp(table_h_margin)
%% figure
% this one works for 2024b

exportgraphics(fig,full_save_path(fig_name),'BackgroundColor','none');
exportgraphics(fig_ACC, full_save_path(fig_name + "_acc"),'BackgroundColor','none','ContentType','vector');
exportgraphics(fig_CCC, full_save_path(fig_name + "_ccc"),'BackgroundColor','none','ContentType','vector');
% saveas(fig, full_save_path(fig_name));
% does not work for png, the margin is still there.
% exportfig_rgb(gcf, char(full_save_path(fig_name)), 'width',18, 'height',9, 'fontmode','fixed', 'fontsize',15, 'color', 'rgb');
% need some additional exe which may is not available.
% eps2pdf(full_save_path(fig_name))
% exportfig_rgb(gcf, char(full_save_path(fig_name)), 'format', 'png', 'width',18, 'height',9, 'fontmode','fixed', 'fontsize',15, 'color', 'rgb');
%% check acceleration smoothing
% figure(2);hold on;grid on; box on;
% plot(CCC_filter_1.profile.time(1:end-1), CCC_filter_1.profile.A1, "b")
% A1_filtered = sgolayfilt(CCC_filter_1.profile.A1, 3, 21);
% plot(CCC_filter_1.profile.time(1:end-1), A1_filtered, "r--")
