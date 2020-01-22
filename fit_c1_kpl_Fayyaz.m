clear all
close all
clc

% read dynamics from MNova text output
%dynamics_MNova = dlmread('C1pyr_20170307.txt');
%[num_fid,num_peak_temp] = size(dynamics_MNova);
%num_peak = (num_peak_temp-1)/4;

dynamics = xlsread("C:\Users\Owner\Downloads\JW355LacPyr");

peak_order = [1 2];

% set kPL estimation parameters
params_est.kPL = 0.02;
params_fix.R1P = 1/30;
params_fix.R1L = 1/25;

pyr_flip =  [4,5,6,7,8, 9,10,12,14,16, 19,21,23,25,28, 31,34,37,41,90];
lac_flip = 90;
flips = zeros(2,length(pyr_flip));
for i = 1:length(pyr_flip)
flips(1,i) = pyr_flip(i)/180*pi;
flips(2,i) = lac_flip/180*pi;
end
spc_params.timepoint_dt = 3;
plot_scale = 1/2;

% fit for kPL
pyr_c2 = dynamics(:,peak_order(1))';
lac_c2 = dynamics(:,peak_order(2))';

pyr_mod = pyr_c2(1:20);
lac_mod = lac_c2(1:20);

S_data = [pyr_mod;lac_mod];

[params_fit, Sfit_lac, ~, ~] = ...
fit_kPL_initial(S_data, spc_params.timepoint_dt, flips, params_fix, params_est, [], 0);

% plot the fit
fprintf('kPL = % .4f(s-1) \n',params_fit.kPL);
Sfit = [pyr_mod;Sfit_lac];
exchange_fitting_1D_plot([], spc_params, S_data, Sfit, plot_scale);




