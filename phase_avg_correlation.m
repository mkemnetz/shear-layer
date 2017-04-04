%%
clearvars;
close all;
clc;

%%
% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed/F2_F14_processed_trimmed.mat');
if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\shear-layer\F3';
elseif isunix
    preamble = '/Volumes/LaCie/MATLAB/Research/Shear Layer/F3';
else
    error('Problem with PC vs Unix determination for loading data');
end

% fname    = 'F2_F14_processed_trimmed.mat';
fname    = 'F2_F14_processed_refined4.mat';
fullname = fullfile(preamble, 'Processed', fname);

load(fullname);

%%
% phase_m     = cat(3, phase_Final{:}); 
% x_Final_m   = cat(3, x_Final{:}); 
% y_Final_m   = cat(3, y_Final{:}); 

%%
% Rescale x and y
x_Final2_m = x_Final_m.*(15e-3);
x_Final2_m = x_Final2_m./(19.7e-3);

y_Final2_m = y_Final_m.*(15e-3);
y_Final2_m = y_Final2_m./(19.7e-3);

%%
% trim4d         = [41; 153; 68; 132];
if(size(phase_m, 2) == 100)
    trim4d         = [1; 100; 1; 100];
elseif(size(phase_m, 2) == 200)
    trim4d         = [1; 200; 1; 200];
else
    error('Trim is incorrect');
end

phase_4d       = trimMat(phase_m, trim4d);
x_4d           = trimMat(x_Final2_m, trim4d);
y_4d           = trimMat(y_Final2_m, trim4d);
[phase_4d, ~ ] = removeTTP( phase_4d , x_4d, y_4d );

%%
phaseShift      = 110;
phaseStart      = 68;
phase_Final_new = phase_4d(:, :, phaseStart:end);

N  = size(phase_Final_new, 3);
N  = size(phase_Final_new, 3) - mod(N, phaseShift);
NN = N/phaseShift;

%%
C      = zeros(2*size(phase_Final_new(:, :, 1), 1)-1, 2*size(phase_Final_new(:, :, 1), 2)-1, N);
C_mean = zeros(2*size(phase_Final_new(:, :, 1), 1)-1, 2*size(phase_Final_new(:, :, 1), 2)-1, phaseShift);
for i = 1:phaseShift
    a = phase_Final_new(:, :, i:phaseShift:N);
    
    for k = 1:size(a, 3)
        C(:, :, k) = xcorr2(a(:, :, k));
    end
    C_mean(:, :, i) = mean(C, 3);
end

%%
figure(); 
surf(C_mean(:, :, 21), 'EdgeColor', 'none'); 
view(2); 
colormap jet;
title('BL Peak');

figure(); 
surf(C_mean(:, :, 45), 'EdgeColor', 'none'); 
view(2); 
colormap jet;
title('BL Trough');

figure(); 
surf(C_mean(:, :, 76), 'EdgeColor', 'none'); 
view(2); 
colormap jet;
title('BL Peak');

figure(); 
surf(C_mean(:, :, 103), 'EdgeColor', 'none'); 
view(2); 
colormap jet;
title('BL Trough');

figure();
plot(C_mean(100, 100:end, 21));
hold on;
plot(C_mean(100, 100:end, 45));
plot(C_mean(100, 100:end, 76));
plot(C_mean(100, 100:end, 103));
title('Streamwise');
legend('peak', 'trough', 'peak', 'trough');

figure();
plot(C_mean(100:end, 100, 21));
hold on;
plot(C_mean(100:end, 100, 45));
plot(C_mean(100:end, 100, 76));
plot(C_mean(100:end, 100, 103));
title('Spanwise');
legend('peak', 'trough', 'peak', 'trough');



