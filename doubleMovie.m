%%
clearvars;
close all;
clc;

%%
% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed/F2_F14_processed_trimmed.mat');
if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\F3';
elseif isunix
    preamble = '/Volumes/LaCie/MATLAB/Research/Shear Layer/F3';
else
    error('Problem with PC vs Unix determination for loading data');
end

fname    = 'F2_F14_processed_trimmed.mat';
fullname = fullfile(preamble, 'Processed', fname);

load(fullname);

%%
phaseShift      = 110;
phaseStart      = 74;
phase_Final_new = phase_Final(phaseStart:end);

N  = size(phase_Final_new, 2);
N  = size(phase_Final_new, 2) - mod(N, phaseShift);
NN = N/phaseShift;

%%
mean_OPD_phi = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);
var_OPD_phi  = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);

for i = 1:phaseShift
    a = phase_Final_new(i:phaseShift:N);
    a = cat(3,a{:}); 
    
    a_mean = mean(a, 3);
    a_var  = var(a,0, 3);
    mean_OPD_phi(:, :, i) = a_mean;
    var_OPD_phi(:, :, i)  = a_var;
    b = a-a_mean;
end

%%
mean_removed_phase = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), size(phase_Final_new{1}, 3));
for k = 1:NN
    a = phase_Final_new;
    a = cat(3,a{:}); 
    c = a(:, :, 1+(k-1)*phaseShift:k*phaseShift);
    
    mean_removed_phase(:, :, 1+(k-1)*phaseShift:k*phaseShift) = ...
        c - mean_OPD_phi;
end

%%
phase_sampled = cat(3,phase_Final_new{:});
phase_sampled = phase_sampled(:, :, 1:N);

x_Final2 = cellfun(@(x) x.*(15e-3),x_Final,'un',0);
x_Final2 = cellfun(@(x) x./(19.7e-3),x_Final2,'un',0);

y_Final2 = cellfun(@(x) x.*(15e-3),y_Final,'un',0);
y_Final2 = cellfun(@(x) x./(19.7e-3),y_Final2,'un',0);
%%
% nframes = N;
nframes = 1500;
phi               = repmat(linspace(0, 4, phaseShift),1,NN);
mean_OPD_phi_copy = repmat(mean_OPD_phi, 1, 1,NN);
var_OPD_phi_copy  = repmat(var_OPD_phi, 1, 1,NN);
Mov(nframes)        = struct('cdata',[],'colormap',[]);

figure(); 
set(gcf,'units','centimeters','position',[0 0 32 24]);
for i = 1:nframes
    f   = 37125; 
    str = sprintf('$\\phi = $ %0.5f, $t = $ %0.5g ',phi(i), i/f); 
    
    subplot(4, 1, 1);
    surf(x_Final2{1}, y_Final2{1}, phase_sampled(:, :, i), 'EdgeColor', 'none'); 
    view(2);
    colormap jet;
    axis tight;
    caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);
    title({'Original'; str}, 'interpreter', 'latex');
    
    subplot(4, 1, 2);
    surf(x_Final2{1}, y_Final2{1}, mean_OPD_phi_copy(:, :, i), 'EdgeColor', 'none'); 
    view(2);
    colormap jet;
    axis tight;
    caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);
    title('Phase-Average', 'interpreter', 'latex'); 
    
    subplot(4, 1, 3);
    surf(x_Final2{1}, y_Final2{1}, var_OPD_phi_copy(:, :, i), 'EdgeColor', 'none');
    view(2);
    colormap jet;
    axis tight;
    caxis([min(min(var_OPD_phi_copy(:, :, 33))) max(max(var_OPD_phi_copy(:, :, 33)))]);
    title({'Variance'; str}, 'interpreter', 'latex'); 
    
    subplot(4, 1, 4);
    surf(x_Final2{1}, y_Final2{1}, mean_removed_phase(:, :, i), 'EdgeColor', 'none'); 
    view(2);
    colormap jet;
    axis tight;
    caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);
    title({'Phase-Average removed'; str}, 'interpreter', 'latex'); 
    
    Mov(i)=getframe(gcf);
end

writeVideoToFile = 1;
if (writeVideoToFile == 1)
    v = VideoWriter('doubleMovie.avi');
    v.FrameRate = 10;
    open(v);
    writeVideo(v,Mov);
    close(v);
end
