clear all
close all
clc

% read dynamics from MNova text output
dynamics_MNova = dlmread('C1pyr_20170307.txt');
[num_fid,num_peak_temp] = size(dynamics_MNova);
num_peak = (num_peak_temp-1)/4;

dynamics = dynamics_MNova(:,(1:num_peak)*4+1);

peak_order = [2 1 3];

% set kPL estimation parameters
params_est.kPL = 0.02;
params_fix.R1P = 1/30;
params_fix.R1L = 1/25;

pyr_flip = 10;
lac_flip = 10;
flips = [pyr_flip/180*pi*ones(1,size(dynamics,1));lac_flip/180*pi*ones(1,size(dynamics,1))];

spc_params.timepoint_dt = 3;
plot_scale = 1/2;

% fit for kPL
pyr_c2 = dynamics(:,peak_order(1))';
lac_c2 = dynamics(:,peak_order(2))';

S_data = [pyr_c2;lac_c2];

[params_fit, Sfit_lac, ~, ~] = ...
fit_kPL(S_data, spc_params.timepoint_dt, flips, params_fix, params_est, [], 0);

% plot the fit
fprintf('kPL = % .4f(s-1) \n',params_fit.kPL);
Sfit = [pyr_c2;Sfit_lac];
exchange_fitting_1D_plot([], spc_params, S_data, Sfit, plot_scale);




