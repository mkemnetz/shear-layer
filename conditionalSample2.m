%%
clear all;
close all;

%%
if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\F3';
elseif isunix
    preamble = '/Volumes/LaCie/MATLAB/Research/Shear Layer/F3';
else
    error('Problem with PC vs Unix determination for loading data');
end

fname    = 'F2_F14_processed.mat';
fullname = fullfile(preamble, 'Processed', fname);

load(fullname);
%%
N    = size(phase_Final, 2);
trim = [4; 32; 3; 69];

trim_a = trim(1);
trim_b = trim(2);
trim_c = trim(3);
trim_d = trim(4);
parfor i = 1:N
    phase_Final{i} = ...
        phase_Final{i}(trim_a:trim_b, trim_c:trim_d);
    x_Final{i}     = ...
        x_Final{i}(trim_a:trim_b, trim_c:trim_d);
    y_Final{i}     = ...
        y_Final{i}(trim_a:trim_b, trim_c:trim_d); 
end

%%
phaseShift      = 110;
% phaseShift      = 55;
% phaseShift      = 80;

% phaseStart      = 68;
phaseStart      = 74;
phase_Final_new = phase_Final(phaseStart:end);

N  = size(phase_Final_new, 2);
N  = size(phase_Final_new, 2) - mod(N, phaseShift);
NN = N/phaseShift;

%%
m = size(phase_Final_new{1}, 1);
n = size(phase_Final_new{1}, 2);
OPDrms = zeros(phaseShift, NN);
phaseStack = zeros(m, n, phaseShift, NN);
for i = 1:phaseShift
    for k = 1:NN
        phaseStack(:, :, i, k) = phase_Final_new{i+phaseShift*(k-1)};
        
        a      = phase_Final_new{i+phaseShift*(k-1)}(14, :);
        b      = a.^2;
        c      = sum(sum(b))./numel(b);
        OPDrms(i, k) = sqrt(c);
    end
end

%%
OPDrms2 = mean(OPDrms, 2);
figure(); plot(linspace(0, 4, phaseShift), OPDrms2);

% %%
% phaseShift      = 110;
% phaseStart      = 68;
% phase_Final_new = phase_Final(phaseStart:end);
% 
% N  = size(phase_Final_new, 2);
% N  = size(phase_Final_new, 2) - mod(N, phaseShift);
% NN = N/phaseShift;
% 
% %%
% OPDrms = zeros(phaseShift, NN);
% for i = 1:phaseShift
%     for k = 1:NN
%         a      = phase_Final_new{i+phaseShift*(k-1)};
%         b      = a.^2;
%         c      = sum(sum(b))./numel(b);
%         OPDrms(i, k) = sqrt(c);
%     end
% end
% 
% OPDrms2 = mean(OPDrms, 2);
% figure(); plot(linspace(0, 4, phaseShift), OPDrms2);
% set(gca,'TickLabelInterpreter','latex')
% set(gca, 'xtickLabel', {'0', '0.5p', 'p', '3p/2', '2p', '5p/2', '3p', '7p/2', '4p'})

%%

meanOPD = zeros(size(phase_Final_new{1}, 2), phaseShift);
meanOPD_2D = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);

for i = 1:phaseShift
    a = phase_Final_new(i:phaseShift:N);
    a = cat(3,a{:}); 
    
    b = a(15, :, :);
    
    meanOPD(:, i) = mean(b, 3);
    meanOPD_2D(:, :, i) = mean(a, 3);
    
end

varOPD = zeros(size(phase_Final_new{1}, 2), phaseShift);
varOPD_2D = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);
varOPD_1D = zeros(1, phaseShift);

mean_OPD_phi = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);

for i = 1:phaseShift
    a = phase_Final_new(i:phaseShift:N);
    a = cat(3,a{:}); 
    
    a_mean = mean(a, 3);
    mean_OPD_phi(:, :, i) = a_mean;
    a = a-a_mean;
    
    b = a(15, :, :);
    c = a(15, 50, :);
    
    varOPD(:, i) = var(b, 0,3);
    varOPD_1D(i) = var(c, 0,3);
    varOPD_2D(:, :, i) = var(a, 0,3);
    
end

%%
load(fullfile(preamble, 'Processed', 'ylocs.mat'));
load(fullfile(preamble, 'Processed', 'pmat.mat'));
load(fullfile(preamble, 'Processed', 'bavgP.mat'));
load(fullfile(preamble, 'synvelmat.mat'));
load(fullfile(preamble, 'yBL.mat'));

fine = 200; 
l = size(pmat, 2); 
w = size(pmat, 1); 
[Xq,Yq] = meshgrid(linspace(1,l, fine), linspace(1,w, fine)); 
[X,Y]   = meshgrid(1:l, 1:w);

bavgP_fine      = repmat(bavgP, 19, 1);
ylocs_fine      = repmat(ylocs', 1, 50);
bavgP_fine      = interp2(X, Y, bavgP_fine, Xq, Yq, 'cubic');
ylocs_fine      = interp2(X, Y, ylocs_fine, Xq, Yq, 'cubic');
pmat_fine       = interp2(X, Y, pmat, Xq, Yq, 'cubic');

l = size(synvelmat, 2); 
w = size(synvelmat, 1); 
[Xq,Yq] = meshgrid(linspace(1,l, fine), linspace(1,w, fine)); 
[X,Y]   = meshgrid(1:l, 1:w);
yBL_fine        = repmat(yBL', 1, 50);
yBL_fine        = interp2(X, Y, yBL_fine, Xq, Yq, 'cubic');
synvelmat_fine  = interp2(X, Y, synvelmat, Xq, Yq, 'cubic');

x_Final2 = cellfun(@(x) x.*(15e-3),x_Final,'un',0);
x_Final2 = cellfun(@(x) x./(19.7e-3),x_Final2,'un',0);

y_Final2 = cellfun(@(x) x.*(15e-3),y_Final,'un',0);
y_Final2 = cellfun(@(x) x./(19.7e-3),y_Final2,'un',0);

%%
% varOPD_phi = squeeze(mean(mean(varOPD, 2), 1));
varOPD_phi = mean(varOPD_1D, 1);
pmat_new = pmat_fine(27, :)+(max(pmat_fine(27, :)-min(pmat_fine(27, :))));
synvelmat_new = synvelmat_fine(157, :)+(max(pmat_fine(157, :)-min(pmat_fine(157, :))));
% varOPD_phi_new = (varOPD_phi./max(varOPD_phi))*5-4;
pmat_new2 = pmat_new./500 +13;
synvelmat_new2 = synvelmat_new;

varOPD_phi_new = (varOPD_phi./max(varOPD_phi));

figure();
% set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(linspace(0, 4, phaseShift), varOPD_phi_new, 'linewidth', 3);
hold on;
% plot(linspace(0, 4, 50), var_vel_norm);
plot(linspace(0, 4, 200), pmat_new2./max(pmat_new2), ':', 'linewidth', 3);
plot(linspace(0, 4, 200), synvelmat_new2./max(synvelmat_new2) - 0.1, 'g--', 'linewidth', 3);
% set(gca,'ytick',[0,pi/2,pi,(3*pi)/2, 2*pi, (5*pi)/2, 3*pi, (7*pi)/2, 4*pi])
% set(gca,'yticklabel',[' 0'; 'pi/2'; 'pi'; '(3*pi)/2'; '2*pi'; '(5*pi)/2'; '3*pi'; '(7*pi)/2'; '4*pi'])

ax = gca;
ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
set(ax, 'FontSize',15);
legend({'Var(OPD)', 'Pressure', 'Velocity'}, 'Interpreter', 'latex', 'location', 'southwest','FontSize',15);
xlabel('Phase (rad.)');

%%
figure();
errorbar(linspace(0, 4, phaseShift),OPDrms2,varOPD_phi)

%%
figure(); 
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);


subplot(2, 2, 1);
surf(x_Final2{1}, y_Final2{1}, meanOPD_2D(:, :, 1), 'EdgeColor', 'none'); 
view(2);
colormap jet;
title({'$\phi = 0$'; '$\overline{OPD}$'}, 'interpreter', 'latex'); 

subplot(2, 2, 3);
surf(x_Final2{1}, y_Final2{1}, varOPD_2D(:, :, 1), 'EdgeColor', 'none'); 
view(2);
colormap jet;
title('Var(OPD)', 'interpreter', 'latex'); 

subplot(2, 2, 2);
% surf(x_Final2{1}, y_Final2{1}, meanOPD_2D(:, :, 83), 'EdgeColor', 'none'); 
surf(x_Final2{1}, y_Final2{1}, meanOPD_2D(:, :, 33), 'EdgeColor', 'none'); 

view(2);
colormap jet;
title({'$\phi = 3\pi$ '; ' $\overline{OPD}$'}, 'interpreter', 'latex'); 

subplot(2, 2, 4);
% surf(x_Final2{1}, y_Final2{1}, varOPD_2D(:, :, 83), 'EdgeColor', 'none'); 
surf(x_Final2{1}, y_Final2{1}, varOPD_2D(:, :, 33), 'EdgeColor', 'none'); 

view(2);
colormap jet;
title('Var(OPD)', 'interpreter', 'latex'); 


figure(); 
set(gcf,'units','centimeters','position',[0 0 1.2*10 8]);


subplot(2, 2, 1);
surf(x_Final2{1}, y_Final2{1}, phase_Final_new{1}, 'EdgeColor', 'none'); 
view(2);
colormap jet;
xlabel('$x/\delta$', 'interpreter', 'latex');
ylabel('$z/\delta$', 'interpreter', 'latex');
title({'$\phi = 0$'; '$OPD$'}, 'interpreter', 'latex'); 
axis tight;
% caxis([min(min(phase_Final_new{83})) max(max(phase_Final_new{83}))]);
caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);


subplot(2, 2, 3);
surf(x_Final2{1}, y_Final2{1}, phase_Final_new{1+4*phaseShift}, 'EdgeColor', 'none'); 
view(2);
xlabel('$x/\delta$', 'interpreter', 'latex');
ylabel('$z/\delta$', 'interpreter', 'latex');
colormap jet;
axis tight;
% caxis([min(min(phase_Final_new{83})) max(max(phase_Final_new{83}))]);
caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);

% title('$OPD$', 'interpreter', 'latex'); 

subplot(2, 2, 2);
% surf(x_Final2{1}, y_Final2{1},  phase_Final_new{83}, 'EdgeColor', 'none'); 
surf(x_Final2{1}, y_Final2{1},  phase_Final_new{33}, 'EdgeColor', 'none'); 

view(2);
xlabel('$x/\delta$', 'interpreter', 'latex');
ylabel('$z/\delta$', 'interpreter', 'latex');
colormap jet;
axis tight;
% caxis([min(min(phase_Final_new{83})) max(max(phase_Final_new{83}))]);
caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);

title({'$\phi = 3\pi$'; '$OPD$'}, 'interpreter', 'latex'); 

subplot(2, 2, 4);
% surf(x_Final2{1}, y_Final2{1}, phase_Final_new{83+4*phaseShift}, 'EdgeColor', 'none'); 
surf(x_Final2{1}, y_Final2{1}, phase_Final_new{33+4*phaseShift}, 'EdgeColor', 'none'); 

view(2);
xlabel('$x/\delta$', 'interpreter', 'latex');
ylabel('$z/\delta$', 'interpreter', 'latex');
colormap jet;
axis tight;
% caxis([min(min(phase_Final_new{83})) max(max(phase_Final_new{83}))]);
caxis([min(min(phase_Final_new{33})) max(max(phase_Final_new{33}))]);


hp4 = get(subplot(2,2,4),'Position');
h = colorbar;
set(h, 'Position', [hp4(1)+hp4(3)+0.03  hp4(2)  0.1  hp4(2)+hp4(3)*2.1]);

% title('$OPD$', 'interpreter', 'latex'); 

% subplot(4, 2, 5);
% surf(x_Final2{1}, y_Final2{1}, phase_Final_new{1}, 'EdgeColor', 'none'); 
% view(2);
% colormap jet;
% 
% subplot(4, 2, 7);
% surf(x_Final2{1+3*phaseShift}, y_Final2{1+phaseShift}, phase_Final_new{1+7*phaseShift}, 'EdgeColor', 'none'); 
% view(2);
% colormap jet;

%%
figure(); 
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);

i = 1; 
subplot(3, 3, 1);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 20; 
subplot(3, 3, 2);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 34; 
subplot(3, 3, 3);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 48; 
subplot(3, 3, 4);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 61; 
subplot(3, 3, 5);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 76; 
subplot(3, 3, 6);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 89; 
subplot(3, 3, 7);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 103; 
subplot(3, 3, 8);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 

i = 110; 
subplot(3, 3, 9);
surf(mean_OPD_phi(:, :, i), 'EdgeColor', 'none'); 
view(2);
colormap jet;
axis tight;
title(aaaa(i),'interpreter', 'latex'); 
