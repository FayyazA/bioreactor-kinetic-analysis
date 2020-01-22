% 1D kPL fitting wrapper
%
% format input data int 4 x timepoints, with rows:
%   pyruvate value
%   lactate value
%   pyruvate flip angles
%   lactate flip angles
%
% include, in separate rows,
%   TR
%   R1P (estimated decay constant for pyruvate)
%   R1L (estimated decay constant for lactate)
%
% Kurhanewicz, Larson, and Vigneron Labs
% 04.2019
%%

clear all
close all

%addpath /home/plarson/matlab/reconstruction

%% data input - reformat depending on data type
filename = 'JW354.mat';
data=load(filename);

% the size and shape of each of these is dependent on the input data
% a=load('JW354.mat') in the matlab command window to interrogate
TR=data.TR;
Nt=size(data.pyruvate,2);
flips=zeros([2,Nt]);
flips(1,:)=data.pyrflip./90;
flips(2,:)=data.lacflip./90;

S=zeros([2,Nt]);
S(1,:)=data.pyruvate;
S(2,:)=data.lactate;

%%

% these three constants are taken from Peder's logic
R1P_est = 1/25; R1L_est = 1/20; % used to be 1/30 and 1/25
kPL_est = .01;
    
% calculate kpl - normal kpl code broken up a bit to handle data shape

% params_fixed, params_est, and params_fit structs are essential to 
% how Peder sets up his code. see in fit_kpl for more info
clear params_fixed params_est
params_fixed.R1P = R1P_est;
params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;

% fit using the input data
% this is not the current version in the hyperpolarized MRI toolbox
% i augmented to handle the change in initial lactate if we begin fitting
% after the bolus
[params_fit,Sfit, ufit, sigma] = fit_kPL(S, TR, flips, params_fixed, params_est,[],0);

% 95% CI on kpl from nlparci
error=(sigma(2)-sigma(1))/2;
Rsquared= 1-sum((reshape(S(2,:),[1 Nt])-Sfit(1,:)).^2)./sum((reshape(S(2,:),[1 Nt])).^2);


%% figure formatting
% these lines look scary, but they are just for formatting the figure
% there is no additional calculation here
% figure
% subplot(121)
% plot(S(1,:))
% hold on
% plot(S(2,:))
% title('Lactate Fit Over Time')
% hold on
% plot(Sfit,'k--')
% xlabel('Timepoints')
% ylabel('a.u.')
% legend('Pyruvate','Lactate','Lactate Fit')
% 
% fprintf('kPL is: %.3f +- %.3f\nR^2 on the lactate fit is: %.3f\n',params_fit.kPL, error,Rsquared)
% xhalf=Nt/2;
% yhalf=max(S(2,:))*0.5;  
% text(xhalf,yhalf+yhalf*0.1,['k_{PL} = ' num2str(params_fit.kPL,'%.3f') '+-' num2str(0.5*(sigma(1,2)-sigma(1,1)),'%.3f')],'fontsize',8);
% text(xhalf,yhalf,['R^{2} = ' num2str(Rsquared,'%.3f')],'fontsize',8);
% fprintf('kPL is: %.3f +- %.3f\nR^2 on the lactate fit is: %.3f\n',params_fit.kPL, error,Rsquared)
% 
% suptitle=annotation('textbox', [0 0.9 1 0.1], 'String', filename, 'EdgeColor', 'none', 'HorizontalAlignment', 'center')
% suptitle.FontSize=16;
% suptitle.FontWeight='bold';

%% optional - for manually removing points before the bolus

% this section is almost identical to the above, but i've manually
% truncated the input data. you can play with this if you want to start 
% fitting at or after the bolus.

S=S(:,4:end); %used to be 4
flips=flips(:,4:end); % used to be 4
Nt=size(flips,2);

R1P_est = 1/30; R1L_est = 1/25;
kPL_est = .01;
    
%calculate kpl - normal kpl code broken up a bit to handle data shape
%really need to better version control this one..
clear params_fixed params_est
params_fixed.R1P = R1P_est;
params_fixed.R1L = R1L_est;
params_est.kPL = kPL_est;

[params_fit,Sfit, ufit, sigma] = fit_kPL(S, TR, flips, params_fixed, params_est,[],0);

% 95% CI on kpl from nlparci
sigma
error=(sigma(2)-sigma(1))/2
Rsquared= 1-sum((reshape(S(2,:),[1 Nt])-Sfit(1,:)).^2)./sum((reshape(S(2,:),[1 Nt])).^2);
fprintf('kPL is: %.3f +- %.3f\nR^2 on the lactate fit is: %.3f\n',params_fit.kPL, error,Rsquared)

% figure formatting
%subplot(122)
plot(S(1,:))
hold on
plot(S(2,:))
title('Lactate Fit Over Time - after Bolus')
hold on
plot(Sfit,'k--')
xlabel('Timepoints')
legend('Pyruvate','Lactate','Lactate Fit')

xhalf=Nt/2;
yhalf=max(S(2,:))*0.5;  
text(xhalf,yhalf+yhalf*0.1,['k_{PL} = ' num2str(params_fit.kPL,'%.3f') '+-' num2str(0.5*(sigma(1,2)-sigma(1,1)),'%.3f')],'fontsize',8);
text(xhalf,yhalf,['R^{2} = ' num2str(Rsquared,'%.3f')],'fontsize',8);
print(gcf,'-dtiff','-r300',strcat(filename,"FINAL",".tif"))