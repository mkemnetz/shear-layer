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
phaseStart      = 66;
phase_Final_new = phase_4d(:, :, phaseStart:end);

N  = size(phase_Final_new, 3);
N  = size(phase_Final_new, 3) - mod(N, phaseShift);
NN = N/phaseShift;

%%
mean_OPD_phi       = zeros(size(phase_Final_new(:, :, 1), 1), size(phase_Final_new(:, :, 1), 2), phaseShift);
var_OPD_phi        = zeros(size(phase_Final_new(:, :, 1), 1), size(phase_Final_new(:, :, 1), 2), phaseShift);


for i = 1:phaseShift
    a = phase_Final_new(:, :, i:phaseShift:N);
%     a = cat(3,a{:}); 
    
    a_mean = mean(a, 3);
    a_var  = var(a,0, 3);
    mean_OPD_phi(:, :, i) = a_mean;
    var_OPD_phi(:, :, i)  = a_var;
    b = a-a_mean;
end
center_var_OPD_phi =  squeeze(var_OPD_phi(100, 100, :))';

%%
mean_removed_phase = zeros(size(phase_Final_new(:, :, 1), 1), size(phase_Final_new(:, :, 1), 2), size(phase_Final_new(:, :, 1), 3));
for k = 1:NN
    a = phase_Final_new;
%     a = cat(3,a{:}); 
    c = a(:, :, 1+(k-1)*phaseShift:k*phaseShift);
    
    mean_removed_phase(:, :, 1+(k-1)*phaseShift:k*phaseShift) = ...
        c - mean_OPD_phi;
end

%%
% phase_sampled = cat(3,phase_Final_new{:});
% phase_sampled = phase_sampled(:, :, 1:N);
phase_sampled = phase_Final_new(:, :, 1:N);

% x_Final2 = cellfun(@(x) x.*(15e-3),x_Final,'un',0);
% x_Final2 = cellfun(@(x) x./(19.7e-3),x_Final2,'un',0);
% 
% y_Final2 = cellfun(@(x) x.*(15e-3),y_Final,'un',0);
% y_Final2 = cellfun(@(x) x./(19.7e-3),y_Final2,'un',0);

%%
phi = repmat(linspace(0, 4, phaseShift),1,NN);

figure(); 
set(gcf,'units','centimeters','position',[0 0 32 24]);
plot(phi(1:phaseShift), center_var_OPD_phi)
% vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'g','r','g','r','g','r','g','r'})
vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'b','r','r','b','b','r','r','b'}, {'SL', 'BL', 'SL', 'BL','SL', 'BL','SL', 'BL'});
xlim([0 4]);
xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
    9/4 5/2 11/4 3 13/4 7/2 15/4 4])
xticklabels({'0' ,'\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4',...
    '3\pi/2', '7\pi/4', '2\pi', '9\pi/4', '5\pi/2', '11\pi/4' ...
    '3\pi', '13\pi/4', '7\pi/2', '15\pi/4', '4\pi'})
title('Temporal Variance at x = 0, y = 0')