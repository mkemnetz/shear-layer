%% Clear workspace 
clear all;
close all;

%% Load in Data
% Note from Piyush:
% trmsmat.mat: I got this using an old version of the code, but the pictures 
% seemed right so it should be good. These are the values of the residual 
% turbulence at each point. It is a 26 x 50 matrix which should link up 
% with everything else 26 x 50. From this, I subtracted the mean of each 
% row to plot Fig 5.15

% Change variable 'preamble' to be a string pointing to the data poarent
% folder.  uigetdir is good for this.
if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\shear-layer\F3\Processed';
elseif isunix
    preamble = '/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed';
else
    error('Problem with PC vs Unix determination for loading data');
end

load(fullfile(preamble, 'trmsmat.mat'));
load(fullfile(preamble, 'yBL.mat'));
load(fullfile(preamble, 'bavgP.mat'));
load(fullfile(preamble, 'Piyush_colormap.mat'));

%% Define Some constants
%  Only U_inf is used in this code but the others are used in other codes
%  so they have been left in.  The actual Mach number Piyush ran at was
%  M ~ 0.625 not M = 0.6.  He does report this in his dissertation.
r     = 0.9;
R     = 286.9;
gamma = 1.4;
M_inf = 0.625;
U_inf = 209.5;

%% Phase average the Urms
%  Note trmsmat isn't actually the turbulence intensity.  From Piyush:  
%  "...it's u_rms then. Sorry, bad nomenclature on my part in using a "t". 
%  Now that you say that I do remember that, but I had used it too many 
%  times in too many codes that I just left it I think since I always
%  knew what it was..."
phaseAverage_U_tilde_rms     = repmat(mean(trmsmat, 2),1,50);

%% Plot the data
figure();
set(gcf,'units','centimeters','position',[0 0 24 8]);

contourf(bavgP./pi, yBL, (trmsmat-phaseAverage_U_tilde_rms)./U_inf);
hold on;
contour(bavgP./pi, yBL, (trmsmat-phaseAverage_U_tilde_rms)./U_inf);

xlim([0 4]);
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})
colormap(cmap);
colorbar('location','eastoutside');
caxis([-0.04 0.04]);
grid on;
grid minor;
title('$(\tilde{U}_{rms} - \langle \tilde{U}_{rms}\rangle_{\phi})/U_0$', 'interpreter', 'latex');
ylabel('$z/\delta$', 'interpreter', 'latex');
xlabel('$x/\delta$', 'interpreter', 'latex');

%%
figure();
set(gcf,'units','centimeters','position',[0 10.5 24 8]);

contourf(bavgP./pi, yBL, (trmsmat)./U_inf);
hold on;
contour(bavgP./pi, yBL, (trmsmat)./U_inf);

xlim([0 4]);
xticks([0 1/2 1 3/2 2 ... 
    5/2 3 7/2 4])
xticklabels({'0', '\pi/2', '\pi',...
    '3\pi/2', '2\pi', '5\pi/2', ...
    '3\pi', '7\pi/2', '4\pi'})

colormap(cmap);
colorbar('location','eastoutside');
caxis([0 0.09]);
grid on;
grid minor;
title('$\tilde{U}_{rms}/U_0$', 'interpreter', 'latex');
ylabel('$z/\delta$', 'interpreter', 'latex');
xlabel('$x/\delta$', 'interpreter', 'latex');
