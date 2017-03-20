%%
clearvars;
close all;
clc;

%%
% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed/F2_F14_processed_refined.mat');
if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\F3';
elseif isunix
    preamble = '/Volumes/LaCie/MATLAB/Research/Shear Layer/F3';
else
    error('Problem with PC vs Unix determination for loading data');
end

fname    = 'F2_F14_processed_refined2.mat';
fullname = fullfile(preamble, 'Processed', fname);

load(fullname);
%%
phaseShift = 110;
phaseStart = 65;
f          = 37125; 

phase_m    = phase_m(:, :, phaseStart:end);

N  = size(phase_m, 3);
N  = size(phase_m, 3) - mod(N, phaseShift);
NN = N/phaseShift;

phase_m    = phase_m(:, :, 1:N);
%%
% Rescale x and y
x_Final2_m = x_Final_m.*(15e-3);
x_Final2_m = x_Final2_m./(19.7e-3);

y_Final2_m = y_Final_m.*(15e-3);
y_Final2_m = y_Final2_m./(19.7e-3);

%% Re-aperture
% Example of trim usage: A(trim(1):trim(2), trim(3):trim(4))
% ~1 delta
trim1d         = [1; 200; 92; 108];
phase_1d       = trimMat(phase_m, trim1d);
x_1d           = trimMat(x_Final2_m, trim1d);
y_1d           = trimMat(y_Final2_m, trim1d);
[phase_1d, ~ ] = removeTTP( phase_1d , x_1d, y_1d );

% ~2 delta
trim2d         = [1; 200; 84; 121];
phase_2d       = trimMat(phase_m, trim2d);
x_2d           = trimMat(x_Final2_m, trim2d);
y_2d           = trimMat(y_Final2_m, trim2d);
[phase_2d, ~ ] = removeTTP( phase_2d , x_2d, y_2d );

% ~3 delta
trim3d         = [1; 200; 76; 124];
phase_3d       = trimMat(phase_m, trim3d);
x_3d           = trimMat(x_Final2_m, trim3d);
y_3d           = trimMat(y_Final2_m, trim3d);
[phase_3d, ~ ] = removeTTP( phase_3d , x_3d, y_3d );

% ~4 delta
trim4d         = [1; 200; 68; 132];
phase_4d       = trimMat(phase_m, trim4d);
x_4d           = trimMat(x_Final2_m, trim4d);
y_4d           = trimMat(y_Final2_m, trim4d);
[phase_4d, ~ ] = removeTTP( phase_4d , x_4d, y_4d );

%% Compute OPDrms
OPDrms_1d = OPDrms_spatial(phase_1d);
OPDrms_2d = OPDrms_spatial(phase_2d);
OPDrms_3d = OPDrms_spatial(phase_3d);
OPDrms_4d = OPDrms_spatial(phase_4d);

%% Phase Average OPDrms
OPDrms_phi_1d = zeros(1, NN);
OPDrms_phi_2d = zeros(1, NN);
OPDrms_phi_3d = zeros(1, NN);
OPDrms_phi_4d = zeros(1, NN);

for i = 1:phaseShift
    a = OPDrms_1d(i:phaseShift:N);
    b = OPDrms_2d(i:phaseShift:N);
    c = OPDrms_3d(i:phaseShift:N);
    d = OPDrms_4d(i:phaseShift:N);
    
    OPDrms_phi_1d(i) = mean(a, 2);
    OPDrms_phi_2d(i) = mean(b, 2);
    OPDrms_phi_3d(i) = mean(c, 2);
    OPDrms_phi_4d(i) = mean(d, 2);
end

%%
phi                = linspace(0, NN*4, NN*phaseShift);
t                  = (1:N).*(1/f);
OPDrms_phi_1d_copy = repmat(OPDrms_phi_1d, 1, NN);
OPDrms_phi_2d_copy = repmat(OPDrms_phi_2d, 1, NN);
OPDrms_phi_3d_copy = repmat(OPDrms_phi_3d, 1, NN);
OPDrms_phi_4d_copy = repmat(OPDrms_phi_4d, 1, NN);

%% Plotting
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(phi, OPDrms_phi_1d_copy);
hold on;
plot(phi, OPDrms_1d);
title('1');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(phi, OPDrms_phi_2d_copy);
hold on;
plot(phi, OPDrms_2d);
title('2');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(phi, OPDrms_phi_3d_copy);
hold on;
plot(phi, OPDrms_3d);
title('3');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(phi, OPDrms_phi_4d_copy);
hold on;
plot(phi, OPDrms_4d);
title('4');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(phi, OPDrms_phi_1d_copy);
hold on;
plot(phi, OPDrms_phi_2d_copy);
plot(phi, OPDrms_phi_3d_copy);
plot(phi, OPDrms_phi_4d_copy);
title('phase average OPDrms');
legend('1d','2d','3d','4d','Location','NorthEastOutside')



figure(); 
set(gcf,'units','centimeters','position',[0 0 32 24]);

subplot(2, 2, 1);
plot(phi(1:331), OPDrms_phi_1d_copy(1:331));
hold on;
plot(phi(1:331), OPDrms_1d(1:331));
axis tight;
title('1');

subplot(2, 2, 2);
plot(phi(1:331), OPDrms_phi_2d_copy(1:331));
hold on;
plot(phi(1:331), OPDrms_2d(1:331));
axis tight;
title('2');

subplot(2, 2, 3);
plot(phi(1:331), OPDrms_phi_3d_copy(1:331));
hold on;
plot(phi(1:331), OPDrms_3d(1:331));
axis tight;
title('3');

subplot(2, 2, 4);
plot(phi(1:331), OPDrms_phi_4d_copy(1:331));
hold on;
plot(phi(1:331), OPDrms_4d(1:331));
axis tight;
title('4');







