clear all; close all; clc;

%% Load images
filename = 'Perfusionimages.mat';
load(filename);

%This shows you the movie of all frames following the more concentrated
%injection
for i =1:1:60
    figure(1); imagesc(squeeze(recon_TF(:,:,i))); colormap gray; caxis([0 600]); pause(0.25)
end

%% Draw ROI for your input function
% The input function is the ventricular blood pool

figure(102); imshow(abs(recon_AIF(:,:,27)),[0 200]); axis image; colormap gray;
[mask_AIF] = roipoly;

z = size(recon_AIF);
for i = 1:z(3)
    figure(101); imshow(abs(recon_AIF(:,:,i)),[0 200]); axis image; colormap gray;
    [B,L] = bwboundaries(mask_AIF,'holes');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    pause(0.1)
end
%% Draw ROI for your tissue function
% Make sure that you draw your region of interest to only include the
% relevant portion of tissue base don AHA segmentation
figure(102); imagesc(abs(recon_TF(:,:,45))); axis image; colormap gray;
[mask_TF] = roipoly;

z = size(recon_TF);
for i = 1:z(3)
    figure(101); imagesc(abs(recon_TF(:,:,i))); axis image; colormap gray;
    [B,L] = bwboundaries(mask_TF,'holes');
    hold on
    for k = 1:length(B)
        boundary = B{k};
        plot(boundary(:,2), boundary(:,1), 'g', 'LineWidth', 2)
    end
    pause(0.1)
end

%% AIF Signal 
clear images image_mat
for i=1:1:z(3)
    Mimage = abs(recon_AIF(:,:,i));
    Mask_image = double(Mimage).*mask_AIF;
    image_mat(i,:) = Mask_image(:);
    images(:,:,i) = Mimage;
end

a = image_mat(1,:); [row idx val] = find(a); %find indices where you actually have data
zps = zeros(1,length(idx));
for i = idx
    data = image_mat(:,i); %all frames of one pixel
end
b = size(image_mat); av_signal = zeros(1,b(1)); std_signal = zeros(1,b(1));
for j = 1:1:b(1)
    vals = image_mat(j,idx); %only average over masked pixels
    av_signal(j) = mean(vals);
    std_signal(j) = std(vals);
end
nt = z(3);
figure; plot(av_signal);
AIFav_signal = av_signal(1,3:nt);
PD_AIFav_signal = mean(av_signal(1,1));

%% TF Signal

clear images image_mat
for i=1:1:z(3)
    Mimage = abs(recon_TF(:,:,i));
    Mask_image = double(Mimage).*mask_TF;
    image_mat(i,:) = Mask_image(:);
    images(:,:,i) = Mimage;
end

a = image_mat(1,:); [row idx val] = find(a); %find indices where you actually have data
zps = zeros(1,length(idx));
for i = idx
    data = image_mat(:,i); %all frames of one pixel
end
b = size(image_mat); av_signal = zeros(1,b(1)); std_signal = zeros(1,b(1));
for j = 1:1:b(1)
    vals = image_mat(j,idx); %only average over masked pixels
    av_signal(j) = mean(vals);
    std_signal(j) = std(vals);
end

figure; plot(av_signal);
TFav_signal = av_signal(1,3:nt);
PD_TFav_signal = mean(av_signal(1,1));

%% Normalized signal
Norm_signal_AIF = AIFav_signal./mean(PD_AIFav_signal);

TF_FA = 12;         % Flip angle for T1 weighted AIF images in degrees
TF_PD_FA = 5;      % Flip angle for Proton density AIF images in degrees

Norm_signal_TF = TFav_signal./mean(PD_TFav_signal);
Norm_factor = sind(TF_FA)/sind(TF_PD_FA);
Norm_signal_TF = Norm_signal_TF./Norm_factor;

%% Signal to Gd

clc; close all;
x0 = 100;
lb = 1;
ub = Inf;
TI_AIF =  5;    % Saturation time for AIF in ms
TI_TF = 90;

for (i = 1:length(Norm_signal_AIF))
    M  = Norm_signal_AIF(i);
    [x,resnorm,residual(i),exitflag]= lsqcurvefit(@SimpleBloch, x0,TI_AIF, M,lb,ub);
    T1(i) = x/1000;
end

R1 = 1./T1;
R1_baseline = mean(1./T1(1:2));
AIFGdconc = (R1 - R1_baseline)./5.3;


for (i = 1:length(Norm_signal_TF))
    M  = Norm_signal_TF(i);
    [x,resnorm,residual(i),exitflag]= lsqcurvefit(@SimpleBloch, x0,TI_TF, M,lb,ub);
    T1_TF(i) = x/1000;
end
R1 = 1./T1_TF;
R1_baseline = mean(R1(1:2));
TFGdconc = (R1 - R1_baseline)./5.2;

%% Fermi Deconvolution

clc; close all;
time = 0:0.8:0.8*length(AIFGdconc)-1;
xdata(:,1) = time(1:41);
xdata(:,2) = AIFGdconc(1:41);
ydata = TFGdconc(1:41)';

x0 = [0.15 1 0.5 0 0.1];
lb = [0 0 0 0 0];
ub = [0.5 Inf Inf Inf Inf];
[rows cols] = size(ydata);

for i = 1:cols
    [x,resnorm,residual,exitflag]= lsqcurvefit(@fermi, x0, xdata, ydata,lb,ub);
    tf_fermi = fermi(x, xdata);
end
tf_fermi = fermi(x, xdata);
figure;
plot(xdata(:,1),ydata,'.',xdata(:,1),tf_fermi);
flow = x(1).*60*0.95;

figure; plot(AIFGdconc);
