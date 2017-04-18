%% 
clear all;
close all;

% trmsmat.mat: I got this using an old version of the code, but the pictures 
% seemed right so it should be good. These are the values of the residual 
% turbulence at each point. It is a 26 x 50 matrix which should link up 
% with everything else 26 x 50. From this, I subtracted the mean of each 
% row to plot Fig 5.15

if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\shear-layer\F3';
elseif isunix
    preamble = '/Volumes/LaCie/MATLAB/Research/Shear Layer/F3';
else
    error('Problem with PC vs Unix determination for loading data');
end

load(fullfile(preamble, 'Processed', 'trmsmat.mat'));
load(fullfile(preamble, 'Processed', 'yBL.mat'));
load(fullfile(preamble, 'Processed', 'bavgP.mat'));
load(fullfile(preamble, 'Processed', 'Piyush_colormap.mat'));
load(fullfile(preamble, 'Processed', 'Ubar.mat'));
load(fullfile(preamble, 'Processed', 'center_temporal_variance2.mat'));
load(fullfile(preamble, 'Processed', 'synvelmat.mat'));

Ue = csvread(fullfile(preamble, 'Processed', 'Ue.txt'));
Ue(:, 1) = Ue(:, 1)./pi;

P0 = csvread(fullfile(preamble, 'Processed', 'P0.txt'));
P0(:, 1) = P0(:, 1)./pi;

[Y,I]=sort(P0(:,1));
P0 = P0(I, :);

Ubar_Mach_full = csvread(fullfile(preamble, 'Processed', 'Ubar.txt'));
%%
r     = 0.9;
R     = 286.9;
gamma = 1.4;
M_inf = 0.625;
U_inf = 209.5;

%%%%%%%%%%%%%%%
% % Meters to microns
% R     = R*((1e6)^2);
% U_inf = U_inf.*1e6;
% Ubar  = Ubar.* 1e6;


cons  = (r*(gamma - 1)*(M_inf^2));
C     = mean(trmsmat, 2);
C_rep = repmat(C,1,50);

a     = U_inf/M_inf;
Ts    = (a^2)./(gamma*R);

%%
figure();
contourf(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
hold on;
contour(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);

colormap(cmap);
colorbar;
caxis([-0.04 0.04]);
grid on;

%%
figure();
contourf(bavgP./pi, yBL, (synvelmat + repmat(Ubar,1,50))./U_inf);
hold on;
contour(bavgP./pi, yBL, (synvelmat + repmat(Ubar,1,50))./U_inf);

colormap(cmap);
colorbar;
caxis([0 1]);
grid on;

figure();
contourf(bavgP./pi, yBL, (repmat(Ubar,1,50))./U_inf);
hold on;
contour(bavgP./pi, yBL, (repmat(Ubar,1,50))./U_inf);

colormap(cmap);
colorbar;
caxis([0 1]);
grid on;


%%
Ubar_full = Ubar_Mach_full(:, 1).*a;
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(Ubar_full./U_inf, Ubar_Mach_full(:, 2));

grid on;
grid minor;
xlabel('$\overline{U}/U_0$', 'interpreter', 'latex');
ylabel('$y/\delta$', 'interpreter', 'latex');
title('Normalized Mean Velocity, $\overline{U}/U_0$', 'interpreter', 'latex');
%%
figure();
set(gcf,'units','centimeters','position',[0 0 24 8]);
contourf(bavgP./pi, yBL, (trmsmat)./U_inf);
hold on;
contour(bavgP./pi, yBL, (trmsmat)./U_inf);
vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'k','k','k','k','k','k','k','k'}, {'SL', 'BL', 'SL', 'BL','SL', 'BL','SL', 'BL'});
xlim([0 4]);
xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
    9/4 5/2 11/4 3 13/4 7/2 15/4 4])
xticklabels({'0' ,'\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4',...
    '3\pi/2', '7\pi/4', '2\pi', '9\pi/4', '5\pi/2', '11\pi/4' ...
    '3\pi', '13\pi/4', '7\pi/2', '15\pi/4', '4\pi'})

colormap(cmap);
colorbar;
caxis([0 0.08]);
grid on;

% openfig('urms.fig');

%%
figure();
plot(yBL, Ubar./U_inf, '.-');
grid on;

%%
[~, ~, ~, rho_ratio, ~] = flowisentropic(gamma, M_inf, 'mach');


% CHECK THESE VALUES
P_total   = convpres(14.25, 'psi', 'Pa');
T_total   = convtemp(75, 'F', 'K');
rho_total = P_total./(R*T_total);
rho_inf   = rho_ratio.*rho_total;

cons2     = cons.*(Ubar./(U_inf^2));
cons2_rep = repmat(cons2, 1, 50);

% cons2     = cons.*((synvelmat + repmat(Ubar,1,50))./(U_inf^2));
% cons2_rep = cons2;

% cons2     = cons.*((synvelmat)./(U_inf^2));
% cons2_rep = cons2;

% cons3     = 1/((1+((gamma-1)/2)*(M_inf^2)*(1-(Ubar./U_inf).^2)).^2);
% cons3_rep = repmat(cons3', 1, 50);
cons3     = 1;
cons3_rep = repmat(cons3', 1, 50);

rho_rms = rho_inf*cons2_rep.*trmsmat.*cons3_rep;

%%
KGD    = 0.000227;
lambda = 0.12;
delta = 19.7e-3;
y     = yBL.*delta;

OPD_rms_squared = 2*delta*(KGD^2)*trapz(y(1:19),(rho_rms(1:19, :).^2)*lambda);


figure()
set(gcf,'units','centimeters','position',[0 0 32 24]);
plot(bavgP./pi, (sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared))).*0.015);
% plot(bavgP./pi, sqrt(OPD_rms_squared));

fig = openfig('center_temporal_variance_408.fig');
figure(fig);
hold on
% plot(bavgP./pi, (sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared))).*0.015);
plot(bavgP./pi, sqrt(OPD_rms_squared).*1e6);

%%
figure();
% set(gcf,'units','centimeters','position',[0 0 24 8]);
set(gcf,'units','centimeters','position',[0 0 1.85*8 12]);

%%%%%%%%%%%%%%%%%%
s1=subplot(3, 1, 1);
contourf(bavgP./pi, yBL, synvelmat./U_inf, 15);
hold on;
contour(bavgP./pi, yBL, synvelmat./U_inf, 15);
p1 = plot (P0(:, 1), (P0(:, 2)*10)+1, '-.g', 'LineWidth', 2);

xlim([0 4]);
% xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
%     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})

colormap(cmap);
hb = colorbar('location','eastoutside');
w = hb.Position;
w(1) = 0.92;
w(2) = 0.735;
w(3) = 0.02;
hb.Position = w;
% caxis([0 1]);
grid on;
grid minor;
title('$(\tilde{u} - \overline{U})/U_0$', 'interpreter', 'latex');
ylabel('$y/\delta$', 'interpreter', 'latex');
hleg = legend([p1], {'$p/P_0$'}, 'Interpreter', 'latex', 'location', 'southeast');
w2 = [0.79, 0.88, 0.1, 0.07];
hleg.Position = w2;
% contourf(bavgP./pi, yBL, (synvelmat + repmat(Ubar,1,50))./U_inf, 15);
% hold on;
% contour(bavgP./pi, yBL, (synvelmat + repmat(Ubar,1,50))./U_inf, 15);
% p1 = plot (P0(:, 1), (P0(:, 2)*10)+1, '-.g', 'LineWidth', 2);
% 
% xlim([0 4]);
% % xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
% %     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
% xticks([0 1/2 1 3/2 2 ... 
%     5/2 3 7/2 4])
% xticklabels({'0', '\pi/2', '\pi',...
%     '3\pi/2', '2\pi', '5\pi/2', ...
%     '3\pi', '7\pi/2', '4\pi'})
% 
% colormap(cmap);
% hb = colorbar('location','eastoutside');
% w = hb.Position;
% w(1) = 0.92;
% w(2) = 0.735;
% w(3) = 0.02;
% hb.Position = w;
% caxis([0 1]);
% grid on;
% grid minor;
% title('$\tilde{u}/U_0$', 'interpreter', 'latex');
% ylabel('$y/\delta$', 'interpreter', 'latex');
% hleg = legend([p1], {'$p/P_0$'}, 'Interpreter', 'latex', 'location', 'southeast');
% w2 = [0.79, 0.88, 0.1, 0.07];
% hleg.Position = w2;


% plot (Ue(:, 1), Ue(:, 2), 'r');
% hold on;
% plot (P0(:, 1), P0(:, 2), '-.b');
% 
% xlim([0 4]);
% % xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
% %     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
% xticks([0 1/2 1 3/2 2 ... 
%     5/2 3 7/2 4])
% xticklabels({'0', '\pi/2', '\pi',...
%     '3\pi/2', '2\pi', '5\pi/2', ...
%     '3\pi', '7\pi/2', '4\pi'})
% ylabel('$y/\delta$', 'interpreter', 'latex');
% hleg = legend({'$u/U_0$', '$p/P_0$'}, 'Interpreter', 'latex', 'location', 'southeast');
% hlegPos = get(hleg,'position');
% w2 = [0.79, 0.76, 0.1, 0.07];
% hleg.Position = w2;

%%%%%%%%%%%%%%%%%%
s2=subplot(3, 1, 2);
contourf(bavgP./pi, yBL, (trmsmat)./U_inf, 10);
hold on;
contour(bavgP./pi, yBL, (trmsmat)./U_inf, 10);
% vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'k','k','k','k','k','k','k','k'}, {'SL', 'BL', 'SL', 'BL','SL', 'BL','SL', 'BL'});
xlim([0 4]);
% xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
%     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})

colormap(cmap);
hb = colorbar('location','eastoutside');
w = hb.Position;
w(1) = 0.92;
w(3) = 0.02;
hb.Position = w;
caxis([0 0.09]);
grid on;
grid minor;
title('$\tilde{u}_{rms}/U_0$', 'interpreter', 'latex');
ylabel('$y/\delta$', 'interpreter', 'latex');

%%%%%%%%%%%%%%%%%%
s3=subplot(3, 1, 3);
plot(phi(1:110), center_var_OPD_phi)
hold on
plot(bavgP./pi, sqrt(OPD_rms_squared).*1e6,'-.', 'LineWidth', 2);
% vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'g','r','g','r','g','r','g','r'})
% vline([0.75 1.60 2.75 3.75],{'r','b','r','b'}, {'Peak', 'Trough', 'Peak', 'Trough'});
xlim([0 4]);
% xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
%     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
% xticklabels({'0' ,'\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4',...
%     '3\pi/2', '7\pi/4', '2\pi', '9\pi/4', '5\pi/2', '11\pi/4' ...
%     '3\pi', '13\pi/4', '7\pi/2', '15\pi/4', '4\pi'})
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})
grid on;
grid minor;
% title('Temporal Variance at x = 0, y = 0')
hleg = legend({'Var(OPD)', 'SRA'}, 'Interpreter', 'latex', 'location', 'northeast');
xlabel('Phase (rad.)', 'interpreter', 'latex');
title('$OPD_{rms}$', 'interpreter', 'latex');

s1Pos = get(s1,'position');
s2Pos = get(s2,'position');
s3Pos = get(s3,'position');


s3Pos(3)   = 0.82; 
s2Pos(3:4) = [s3Pos(3:4)];
s1Pos(3:4) = [s3Pos(3:4)];

s1Pos(1)   = 0.095;
s1Pos(2)   = s1Pos(2)+0.025; 
s2Pos(1)   = 0.095; 
s3Pos(1)   = 0.095;  
s3Pos(2)   = s3Pos(2)-0.025; 
set(s1,'position',s1Pos);
set(s2,'position',s2Pos);
set(s3,'position',s3Pos);

hlegPos = get(hleg,'position');
w2 = [0.77, 0.23, 0.1, 0.07];

hleg.Position = w2;

% %%
% figure(); 
% set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2, 1, 1);
% contourf(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% hold on;
% contour(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% % ax = gca;
% % ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
% 
% colormap(cmap);
% colorbar;
% caxis([-0.04 0.04]);
% grid on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2, 1, 2);
% plot(bavgP./pi, sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared)));
% ax = gca;
% % ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
% % legend({'Var(OPD)', 'Pressure', 'Velocity'}, 'Interpreter', 'latex', 'location', 'bestoutside');
% xlabel('Phase (rad.)');

%%
% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/plottingData.mat')
% 
% figure(); 
% set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2, 1, 1);
% contourf(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% hold on;
% contour(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% % ax = gca;
% % ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
% 
% colormap(cmap);
% colorbar;
% caxis([-0.04 0.04]);
% grid on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2, 1, 2);
% plot(linspace(0, 4, 110), plotted_varOPD_phi./max(plotted_varOPD_phi));
% hold on;
% plot(linspace(0, 4, 200), plotted_pmat, ':');
% % plot(linspace(0, 4, 200), plotted_synvelmat, 'g--');
% % plot(bavgP./pi, ((sqrt(OPD)./max(sqrt(OPD)))+0.6)./max(((sqrt(OPD)./max(sqrt(OPD)))+0.6)));
% plot(bavgP./pi, sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared)));
% ax = gca;
% % ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
% legend({'Var(OPD)', 'Pressure', 'Velocity'}, 'Interpreter', 'latex', 'location', 'bestoutside');
% xlabel('Phase (rad.)');

%%

load(fullfile(preamble, 'Processed', 'SRA_utilde.mat'));
figure();
% set(gcf,'units','centimeters','position',[0 0 24 8]);
set(gcf,'units','centimeters','position',[0 0 1.85*8 12]);

%%%%%%%%%%%%%%%%%%
s1=subplot(3, 1, 1);
contourf(bavgP./pi, yBL, synvelmat./U_inf, 15);
hold on;
contour(bavgP./pi, yBL, synvelmat./U_inf, 15);
p1 = plot (P0(:, 1), (P0(:, 2)*10)+1, '-.g', 'LineWidth', 2);

xlim([0 4]);
% xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
%     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})

colormap(cmap);
hb = colorbar('location','eastoutside');
w = hb.Position;
w(1) = 0.92;
w(2) = 0.735;
w(3) = 0.02;
hb.Position = w;
% caxis([0 1]);
grid on;
grid minor;
title('$(\tilde{u} - \overline{U})/U_0$', 'interpreter', 'latex');
ylabel('$y/\delta$', 'interpreter', 'latex');
hleg = legend([p1], {'$p/P_0$'}, 'Interpreter', 'latex', 'location', 'southeast');
w2 = [0.79, 0.88, 0.1, 0.07];
hleg.Position = w2;
% contourf(bavgP./pi, yBL, (synvelmat + repmat(Ubar,1,50))./U_inf, 15);
% hold on;
% contour(bavgP./pi, yBL, (synvelmat + repmat(Ubar,1,50))./U_inf, 15);
% p1 = plot (P0(:, 1), (P0(:, 2)*10)+1, '-.g', 'LineWidth', 2);
% 
% xlim([0 4]);
% % xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
% %     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
% xticks([0 1/2 1 3/2 2 ... 
%     5/2 3 7/2 4])
% xticklabels({'0', '\pi/2', '\pi',...
%     '3\pi/2', '2\pi', '5\pi/2', ...
%     '3\pi', '7\pi/2', '4\pi'})
% 
% colormap(cmap);
% hb = colorbar('location','eastoutside');
% w = hb.Position;
% w(1) = 0.92;
% w(2) = 0.735;
% w(3) = 0.02;
% hb.Position = w;
% caxis([0 1]);
% grid on;
% grid minor;
% title('$\tilde{u}/U_0$', 'interpreter', 'latex');
% ylabel('$y/\delta$', 'interpreter', 'latex');
% hleg = legend([p1], {'$p/P_0$'}, 'Interpreter', 'latex', 'location', 'southeast');
% w2 = [0.79, 0.88, 0.1, 0.07];
% hleg.Position = w2;


% plot (Ue(:, 1), Ue(:, 2), 'r');
% hold on;
% plot (P0(:, 1), P0(:, 2), '-.b');
% 
% xlim([0 4]);
% % xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
% %     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
% xticks([0 1/2 1 3/2 2 ... 
%     5/2 3 7/2 4])
% xticklabels({'0', '\pi/2', '\pi',...
%     '3\pi/2', '2\pi', '5\pi/2', ...
%     '3\pi', '7\pi/2', '4\pi'})
% ylabel('$y/\delta$', 'interpreter', 'latex');
% hleg = legend({'$u/U_0$', '$p/P_0$'}, 'Interpreter', 'latex', 'location', 'southeast');
% hlegPos = get(hleg,'position');
% w2 = [0.79, 0.76, 0.1, 0.07];
% hleg.Position = w2;

%%%%%%%%%%%%%%%%%%
s2=subplot(3, 1, 2);
contourf(bavgP./pi, yBL, (trmsmat)./U_inf, 10);
hold on;
contour(bavgP./pi, yBL, (trmsmat)./U_inf, 10);
% vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'k','k','k','k','k','k','k','k'}, {'SL', 'BL', 'SL', 'BL','SL', 'BL','SL', 'BL'});
xlim([0 4]);
% xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
%     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})

colormap(cmap);
hb = colorbar('location','eastoutside');
w = hb.Position;
w(1) = 0.92;
w(3) = 0.02;
hb.Position = w;
caxis([0 0.09]);
grid on;
grid minor;
title('$\tilde{u}_{rms}/U_0$', 'interpreter', 'latex');
ylabel('$y/\delta$', 'interpreter', 'latex');

%%%%%%%%%%%%%%%%%%
s3=subplot(3, 1, 3);
plot(phi(1:110), center_var_OPD_phi)
hold on
plot(bavgP./pi, sqrt(OPD_rms_squared).*1e6,'-.', 'LineWidth', 2);
plot(bavgP./pi, SRA_utilde,'-.', 'LineWidth', 2);
% vline([0.5 0.75 1.35 1.60 2.4 2.75 3.35 3.75],{'g','r','g','r','g','r','g','r'})
% vline([0.75 1.60 2.75 3.75],{'r','b','r','b'}, {'Peak', 'Trough', 'Peak', 'Trough'});
xlim([0 4]);
% xticks([0 1/4 1/2 3/4 1 5/4 3/2 7/4 2 ... 
%     9/4 5/2 11/4 3 13/4 7/2 15/4 4])
% xticklabels({'0' ,'\pi/4', '\pi/2', '3\pi/4', '\pi', '5\pi/4',...
%     '3\pi/2', '7\pi/4', '2\pi', '9\pi/4', '5\pi/2', '11\pi/4' ...
%     '3\pi', '13\pi/4', '7\pi/2', '15\pi/4', '4\pi'})
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})
grid on;
grid minor;
% title('Temporal Variance at x = 0, y = 0')
hleg = legend({'Var(OPD)', 'SRA'}, 'Interpreter', 'latex', 'location', 'northeast');
xlabel('Phase (rad.)', 'interpreter', 'latex');
title('$OPD_{rms}$', 'interpreter', 'latex');

s1Pos = get(s1,'position');
s2Pos = get(s2,'position');
s3Pos = get(s3,'position');


s3Pos(3)   = 0.82; 
s2Pos(3:4) = [s3Pos(3:4)];
s1Pos(3:4) = [s3Pos(3:4)];

s1Pos(1)   = 0.095;
s1Pos(2)   = s1Pos(2)+0.025; 
s2Pos(1)   = 0.095; 
s3Pos(1)   = 0.095;  
s3Pos(2)   = s3Pos(2)-0.025; 
set(s1,'position',s1Pos);
set(s2,'position',s2Pos);
set(s3,'position',s3Pos);

hlegPos = get(hleg,'position');
w2 = [0.77, 0.23, 0.1, 0.07];

hleg.Position = w2;

% %%
% figure(); 
% set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2, 1, 1);
% contourf(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% hold on;
% contour(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% % ax = gca;
% % ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
% 
% colormap(cmap);
% colorbar;
% caxis([-0.04 0.04]);
% grid on;
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(2, 1, 2);
% plot(bavgP./pi, sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared)));
% ax = gca;
% % ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
% % legend({'Var(OPD)', 'Pressure', 'Velocity'}, 'Interpreter', 'latex', 'location', 'bestoutside');
% xlabel('Phase (rad.)');
