%% 
clear all;
close all;

% trmsmat.mat: I got this using an old version of the code, but the pictures 
% seemed right so it should be good. These are the values of the residual 
% turbulence at each point. It is a 26 x 50 matrix which should link up 
% with everything else 26 x 50. From this, I subtracted the mean of each 
% row to plot Fig 5.15



load('/Volumes/LaCie/MATLAB/Research/Shear Layer/trmsmat.mat');

load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/yBL.mat');
load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed/bavgP.mat');
load('/Volumes/LaCie/MATLAB/Research/Shear Layer/Piyush_colormap.mat');
load('/Volumes/LaCie/MATLAB/Research/Shear Layer/Ubar.mat');
%%
r     = 0.9;
R     = 286.9;
gamma = 1.4;
M_inf = 0.6;
U_inf = 209.5;


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

openfig('urms.fig');

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

rho_rms = rho_inf*cons2_rep.*trmsmat;

%%
KGD    = 0.000227;
lambda = 0.1;

OPD_rms_squared = 2*(KGD^2)*trapz(yBL(1:21),(rho_rms(1:21, :).^2)*lambda);
figure()
plot(bavgP./pi, sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared)));

%%
load('/Volumes/LaCie/MATLAB/Research/Shear Layer/plottingData.mat')

figure(); 
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);

%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 1, 1);
contourf(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
hold on;
contour(bavgP./pi, yBL, (trmsmat-C_rep)./U_inf);
% ax = gca;
% ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';

colormap(cmap);
colorbar;
caxis([-0.04 0.04]);
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(2, 1, 2);
plot(linspace(0, 4, 110), plotted_varOPD_phi./max(plotted_varOPD_phi));
hold on;
plot(linspace(0, 4, 200), plotted_pmat, ':');
% plot(linspace(0, 4, 200), plotted_synvelmat, 'g--');
% plot(bavgP./pi, ((sqrt(OPD)./max(sqrt(OPD)))+0.6)./max(((sqrt(OPD)./max(sqrt(OPD)))+0.6)));
plot(bavgP./pi, sqrt(OPD_rms_squared)./max(sqrt(OPD_rms_squared)));
ax = gca;
% ax.XAxis.TickLabels = {' 0'; '\pi/2'; '\pi'; '3\pi/2'; '2\pi'; '5\pi/2'; '3\pi'; '7\pi/2'; '4\pi'}';
legend({'Var(OPD)', 'Pressure', 'Velocity'}, 'Interpreter', 'latex', 'location', 'bestoutside');
xlabel('Phase (rad.)');
