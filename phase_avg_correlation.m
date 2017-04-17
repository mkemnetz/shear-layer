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


%% Re-aperture
% Example of trim usage: A(trim(1):trim(2), trim(3):trim(4))


% % trim4d         = [41; 153; 68; 132];
% if(size(phase_m, 2) == 100)
%     trim4d         = [1; 100; 1; 100];
% elseif(size(phase_m, 2) == 200)
%     trim4d         = [1; 200; 1; 200];
% else
%     error('Trim is incorrect');
% end

if(size(phase_m, 2) == 100)
    trim2d         = [1; 100; 42; 58];
    trim3d         = [1; 100; 38; 62];
    trim4d         = [1; 100; 34; 66];
elseif(size(phase_m, 2) == 200)
    trim2d         = [1; 200; 1; 200];
    trim3d         = [1; 200; 1; 200];
    trim4d         = [1; 200; 1; 200];
    warning('Re-aperture not implemented for 200x200 data.  Using Full  2D wavefront')
else
    error('Trim is incorrect');
end
% ~2 delta
phase_2d       = trimMat(phase_m, trim2d);
x_2d           = trimMat(x_Final2_m, trim2d);
y_2d           = trimMat(y_Final2_m, trim2d);
[phase_2d, ~ ] = removeTTP( phase_2d , x_2d, y_2d );

% ~3 delta
phase_3d       = trimMat(phase_m, trim3d);
x_3d           = trimMat(x_Final2_m, trim3d);
y_3d           = trimMat(y_Final2_m, trim3d);
[phase_3d, ~ ] = removeTTP( phase_3d , x_3d, y_3d );

% ~4 delta
phase_4d       = trimMat(phase_m, trim4d);
x_4d           = trimMat(x_Final2_m, trim4d);
y_4d           = trimMat(y_Final2_m, trim4d);
[phase_4d, ~ ] = removeTTP( phase_4d , x_4d, y_4d );

%%
phaseShift      = 110;
phaseStart      = 68;

phase_Final_new_2d = phase_2d(:, :, phaseStart:end);
phase_Final_new_3d = phase_3d(:, :, phaseStart:end);
phase_Final_new_4d = phase_4d(:, :, phaseStart:end);

N  = size(phase_Final_new_4d, 3);
N  = size(phase_Final_new_4d, 3) - mod(N, phaseShift);
NN = N/phaseShift;

%%
C_2d      = zeros(2*size(phase_Final_new_2d(:, :, 1), 1)-1, 2*size(phase_Final_new_2d(:, :, 1), 2)-1, N);
C_mean_2d = zeros(2*size(phase_Final_new_2d(:, :, 1), 1)-1, 2*size(phase_Final_new_2d(:, :, 1), 2)-1, phaseShift);

C_3d      = zeros(2*size(phase_Final_new_3d(:, :, 1), 1)-1, 2*size(phase_Final_new_3d(:, :, 1), 2)-1, N);
C_mean_3d = zeros(2*size(phase_Final_new_3d(:, :, 1), 1)-1, 2*size(phase_Final_new_3d(:, :, 1), 2)-1, phaseShift);

C_4d      = zeros(2*size(phase_Final_new_4d(:, :, 1), 1)-1, 2*size(phase_Final_new_4d(:, :, 1), 2)-1, N);
C_mean_4d = zeros(2*size(phase_Final_new_4d(:, :, 1), 1)-1, 2*size(phase_Final_new_4d(:, :, 1), 2)-1, phaseShift);

h     = waitbar(0,'Please wait...');
steps = phaseShift;
for i = 1:phaseShift
    waitbar(i / steps);
    a_2d = phase_Final_new_2d(:, :, i:phaseShift:N);
    a_3d = phase_Final_new_3d(:, :, i:phaseShift:N);
    a_4d = phase_Final_new_4d(:, :, i:phaseShift:N);
    
    for k = 1:size(a_2d, 3)
        C_2d(:, :, k) = xcorr2(a_2d(:, :, k), a_2d(:, :, k));
    end
    C_mean_2d(:, :, i) = mean(C_2d, 3);
    
    for k = 1:size(a_3d, 3)
        C_3d(:, :, k) = xcorr2(a_3d(:, :, k), a_3d(:, :, k));
    end
    C_mean_3d(:, :, i) = mean(C_3d, 3);
    
    for k = 1:size(a_4d, 3)
        C_4d(:, :, k) = xcorr2(a_4d(:, :, k), a_4d(:, :, k));
    end
    C_mean_4d(:, :, i) = mean(C_4d, 3);
end
close(h) 

% 2d
C_complete_2d      = zeros(2*size(phase_Final_new_2d(:, :, 1), 1)-1, 2*size(phase_Final_new_2d(:, :, 1), 2)-1, N);
for k = 1:N
    C_complete_2d(:, :, k) = xcorr2(phase_Final_new_2d(:, :, k));
end
C_mean_complete_2d = mean(C_complete_2d, 3);

% 3d
C_complete_3d      = zeros(2*size(phase_Final_new_3d(:, :, 1), 1)-1, 2*size(phase_Final_new_3d(:, :, 1), 2)-1, N);
for k = 1:N
    C_complete_3d(:, :, k) = xcorr2(phase_Final_new_3d(:, :, k));
end
C_mean_complete_3d = mean(C_complete_3d, 3);

% 4d
C_complete_4d      = zeros(2*size(phase_Final_new_4d(:, :, 1), 1)-1, 2*size(phase_Final_new_4d(:, :, 1), 2)-1, N);
for k = 1:N
    C_complete_4d(:, :, k) = xcorr2(phase_Final_new_4d(:, :, k));
end
C_mean_complete_4d = mean(C_complete_4d, 3);

%%

% figure(); 
% surf(C_mean(:, :, 21)./C_mean(100, 100, 21), 'EdgeColor', 'none'); 
% view(2); 
% colormap jet;
% title('BL Peak');
% 
% figure(); 
% surf(C_mean(:, :, 45)./C_mean(100, 100, 45), 'EdgeColor', 'none'); 
% view(2); 
% colormap jet;
% title('BL Trough');
% 
% figure(); 
% surf(C_mean(:, :, 76)./C_mean(100, 100, 76), 'EdgeColor', 'none'); 
% view(2); 
% colormap jet;
% title('BL Peak');
% 
% figure(); 
% surf(C_mean(:, :, 103)./C_mean(100, 100, 103), 'EdgeColor', 'none'); 
% view(2); 
% colormap jet;
% title('BL Trough');
% 
% figure();
% plot(C_mean(100, 100:end, 21)./C_mean(100, 100, 21));
% hold on;
% plot(C_mean(100, 100:end, 45)./C_mean(100, 100, 45));
% plot(C_mean(100, 100:end, 76)./C_mean(100, 100, 76));
% plot(C_mean(100, 100:end, 103)./C_mean(100, 100, 103));
% title('Streamwise');
% legend('peak', 'trough', 'peak', 'trough');
% 
% figure();
% plot(C_mean(100:end, 100, 21)./C_mean(100, 100, 21));
% hold on;
% plot(C_mean(100:end, 100, 45)./C_mean(100, 100, 45));
% plot(C_mean(100:end, 100, 76)./C_mean(100, 100, 76));
% plot(C_mean(100:end, 100, 103)./C_mean(100, 100, 103));
% title('Spanwise');
% legend('peak', 'trough', 'peak', 'trough');

%%
% 2d
C_peak_st_2d          = mean([C_mean_2d(100, 17:end, 21)./C_mean_2d(100, 17, 21); C_mean_2d(100, 17:end, 76)./C_mean_2d(100, 17, 76)], 1);
C_trough_st_2d        = mean([C_mean_2d(100, 17:end, 45)./C_mean_2d(100, 17, 45); C_mean_2d(100, 17:end, 103)./C_mean_2d(100, 17, 103)], 1);

C_peak_sp_2d          = mean([C_mean_2d(100:end, 17, 21)'./C_mean_2d(100, 17, 21)'; C_mean_2d(100:end, 17, 76)'./C_mean_2d(100, 17, 76)'], 1);
C_trough_sp_2d        = mean([C_mean_2d(100:end, 17, 45)'./C_mean_2d(100, 17, 45)'; C_mean_2d(100:end, 17, 103)'./C_mean_2d(100, 17, 103)'], 1);
x_2d_shift            = x_2d + -1*x_2d(1, 1, 1);
y_2d_shift            = y_2d + -1*y_2d(1, 1, 1);

C_mean_complete_st_2d = C_mean_complete_2d(100, 17:end)./C_mean_complete_2d(100, 17);
C_mean_complete_sp_2d = C_mean_complete_2d(100:end, 17)'./C_mean_complete_2d(100, 17)';

% 3d
C_peak_st_3d          = mean([C_mean_3d(100, 25:end, 21)./C_mean_3d(100, 25, 21); C_mean_3d(100, 25:end, 76)./C_mean_3d(100, 25, 76)], 1);
C_trough_st_3d        = mean([C_mean_3d(100, 25:end, 45)./C_mean_3d(100, 25, 45); C_mean_3d(100, 25:end, 103)./C_mean_3d(100, 25, 103)], 1);

C_peak_sp_3d          = mean([C_mean_3d(100:end, 25, 21)'./C_mean_3d(100, 25, 21)'; C_mean_3d(100:end, 25, 76)'./C_mean_3d(100, 25, 76)'], 1);
C_trough_sp_3d        = mean([C_mean_3d(100:end, 25, 45)'./C_mean_3d(100, 25, 45)'; C_mean_3d(100:end, 25, 103)'./C_mean_3d(100, 25, 103)'], 1);
x_3d_shift            = x_3d + -1*x_3d(1, 1, 1);
y_3d_shift            = y_3d + -1*y_3d(1, 1, 1);

C_mean_complete_st_3d = C_mean_complete_3d(100, 25:end)./C_mean_complete_3d(100, 25);
C_mean_complete_sp_3d = C_mean_complete_3d(100:end, 25)'./C_mean_complete_3d(100, 25)';

% 4d
C_peak_st_4d          = mean([C_mean_4d(100, 33:end, 21)./C_mean_4d(100, 33, 21); C_mean_4d(100, 33:end, 76)./C_mean_4d(100, 33, 76)], 1);
C_trough_st_4d        = mean([C_mean_4d(100, 33:end, 45)./C_mean_4d(100, 33, 45); C_mean_4d(100, 33:end, 103)./C_mean_4d(100, 33, 103)], 1);

C_peak_sp_4d          = mean([C_mean_4d(100:end, 33, 21)'./C_mean_4d(100, 33, 21)'; C_mean_4d(100:end, 33, 76)'./C_mean_4d(100, 33, 76)'], 1);
C_trough_sp_4d        = mean([C_mean_4d(100:end, 33, 45)'./C_mean_4d(100, 33, 45)'; C_mean_4d(100:end, 33, 103)'./C_mean_4d(100, 33, 103)'], 1);
x_4d_shift            = x_4d + -1*x_4d(1, 1, 1);
y_4d_shift            = y_4d + -1*y_4d(1, 1, 1);

C_mean_complete_st_4d = C_mean_complete_4d(100, 33:end)./C_mean_complete_4d(100, 33);
C_mean_complete_sp_4d = C_mean_complete_4d(100:end, 33)'./C_mean_complete_4d(100, 33)';

%%
% 2d
%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(x_2d_shift(1, :, 1), C_peak_st_2d, '-.r');
hold on;
h1=plot(x_2d_shift(1, 1:3:end, 1), C_peak_st_2d(1:3:end), 'rd');

plot(x_2d_shift(1, :, 1), C_trough_st_2d, '--b');
h2=plot(x_2d_shift(1, 1:3:end, 1), C_trough_st_2d(1:3:end), 'b^');
plot(x_2d_shift(1, :, 1), C_mean_complete_st_2d, '-g');
h3=plot(x_2d_shift(1, 1:3:end, 1), C_mean_complete_st_2d(1:3:end), 'gx');
title('Streamwise, $2 \delta$', 'interpreter', 'latex');
hleg=legend([h1 h2 h3],{'peak','trough', 'average'});
set(hleg, 'interpreter', 'latex');
xlim([0 2]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_x(\Delta x, \Delta z = 0)$', 'interpreter', 'latex');
grid on;
grid minor;

%%%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(y_2d_shift(:, 1, 1), C_peak_sp_2d, '-.r');
hold on;
h1=plot(y_2d_shift(1:3:end, 1, 1), C_peak_sp_2d(1:3:end), 'rd');

plot(y_2d_shift(:, 1, 1), C_trough_sp_2d, '--b');
h2=plot(y_2d_shift(1:3:end, 1, 1), C_trough_sp_2d(1:3:end), 'b^');
plot(y_2d_shift(:, 1, 1), C_mean_complete_sp_2d, '-g');
h3=plot(y_2d_shift(1:3:end, 1, 1), C_mean_complete_sp_2d(1:3:end), 'gx');
title('Spanwise, $2 \delta$', 'interpreter', 'latex');
hleg=legend([h1 h2 h3],{'peak','trough', 'average'});
set(hleg, 'interpreter', 'latex', 'interpreter', 'latex');
xlim([0 5]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_z(\Delta x = 0, \Delta z)$', 'interpreter', 'latex');
grid on;
grid minor;

%%
% 3d
%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(x_3d_shift(1, :, 1), C_peak_st_3d, '-.r');
hold on;
h1=plot(x_3d_shift(1, 1:3:end, 1), C_peak_st_3d(1:3:end), 'rd');

plot(x_3d_shift(1, :, 1), C_trough_st_3d, '--b');
h2=plot(x_3d_shift(1, 1:3:end, 1), C_trough_st_3d(1:3:end), 'b^');
plot(x_3d_shift(1, :, 1), C_mean_complete_st_3d, '-g');
h3=plot(x_3d_shift(1, 1:3:end, 1), C_mean_complete_st_3d(1:3:end), 'gx');
title('Streamwise, $3 \delta$', 'interpreter', 'latex');
hleg=legend([h1 h2 h3],{'peak','trough', 'average'});
set(hleg, 'interpreter', 'latex');
xlim([0 3]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_x(\Delta x, \Delta z = 0)$', 'interpreter', 'latex');
grid on;
grid minor;

%%%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(y_3d_shift(:, 1, 1), C_peak_sp_3d, '-.r');
hold on;
h1=plot(y_3d_shift(1:3:end, 1, 1), C_peak_sp_3d(1:3:end), 'rd');

plot(y_3d_shift(:, 1, 1), C_trough_sp_3d, '--b');
h2=plot(y_3d_shift(1:3:end, 1, 1), C_trough_sp_3d(1:3:end), 'b^');
plot(y_3d_shift(:, 1, 1), C_mean_complete_sp_3d, '-g');
h3=plot(y_3d_shift(1:3:end, 1, 1), C_mean_complete_sp_3d(1:3:end), 'gx');
title('Spanwise, $3 \delta$', 'interpreter', 'latex');
hleg=legend([h1 h2 h3],{'peak','trough', 'average'});
set(hleg, 'interpreter', 'latex', 'interpreter', 'latex');
xlim([0 5]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_z(\Delta x = 0, \Delta z)$', 'interpreter', 'latex');
grid on;
grid minor;

%%
% 4d
%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(x_4d_shift(1, :, 1), C_peak_st_4d, '-.r');
hold on;
h1=plot(x_4d_shift(1, 1:3:end, 1), C_peak_st_4d(1:3:end), 'rd');

plot(x_4d_shift(1, :, 1), C_trough_st_4d, '--b');
h2=plot(x_4d_shift(1, 1:3:end, 1), C_trough_st_4d(1:3:end), 'b^');
plot(x_4d_shift(1, :, 1), C_mean_complete_st_4d, '-g');
h3=plot(x_4d_shift(1, 1:3:end, 1), C_mean_complete_st_4d(1:3:end), 'gx');
title('Streamwise, $Ap/\delta = 4$', 'interpreter', 'latex');
hleg=legend([h1 h2 h3],{'Local max','Local min', 'average'});
set(hleg, 'interpreter', 'latex');
xlim([0 4]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_x(\Delta x, \Delta z = 0)$', 'interpreter', 'latex');
grid on;
grid minor;

%%%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(y_4d_shift(:, 1, 1), C_peak_sp_4d, '-.r');
hold on;
h1=plot(y_4d_shift(1:3:end, 1, 1), C_peak_sp_4d(1:3:end), 'rd');

plot(y_4d_shift(:, 1, 1), C_trough_sp_4d, '--b');
h2=plot(y_4d_shift(1:3:end, 1, 1), C_trough_sp_4d(1:3:end), 'b^');
plot(y_4d_shift(:, 1, 1), C_mean_complete_sp_4d, '-g');
h3=plot(y_4d_shift(1:3:end, 1, 1), C_mean_complete_sp_4d(1:3:end), 'gx');
title('Spanwise, $Ap/\delta = 4$', 'interpreter', 'latex');
hleg=legend([h1 h2 h3],{'Local max','Local min', 'Average'});
set(hleg, 'interpreter', 'latex', 'interpreter', 'latex');
xlim([0 4]);
xlabel('$\Delta y/\delta$', 'interpreter', 'latex');
ylabel('$\rho_z(\Delta x = 0, \Delta z)$', 'interpreter', 'latex');
grid on;
grid minor;

%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(x_2d_shift(1, :, 1), C_peak_st_2d, '-.r');
hold on;
h1=plot(x_2d_shift(1, 1:3:end, 1), C_peak_st_2d(1:3:end), 'rd');

plot(x_2d_shift(1, :, 1), C_trough_st_2d, '--b');
h2=plot(x_2d_shift(1, 1:3:end, 1), C_trough_st_2d(1:3:end), 'b^');

plot(x_3d_shift(1, :, 1), C_peak_st_3d, '-.g');
hold on;
h3=plot(x_3d_shift(1, 1:3:end, 1), C_peak_st_3d(1:3:end), 'gp');

plot(x_3d_shift(1, :, 1), C_trough_st_3d, '--k');
h4=plot(x_3d_shift(1, 1:3:end, 1), C_trough_st_3d(1:3:end), 'ks');

title('Streamwise', 'interpreter', 'latex');
hleg=legend([h1 h2 h3 h4],{'peak, $Ap/\delta = 2$','trough, $Ap/\delta = 2$', 'peak, $Ap/\delta = 3$','trough, $Ap/\delta = 3$'});
set(hleg, 'interpreter', 'latex');
xlim([0 3]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_x(\Delta x, \Delta z = 0)$', 'interpreter', 'latex');
grid on;
grid minor;

%%%%%%%%%%%%%%%%%%%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
plot(y_2d_shift(:, 1, 1), C_peak_sp_2d, '-.r');
hold on;
h1=plot(y_2d_shift(1:3:end, 1, 1), C_peak_sp_2d(1:3:end), 'rd');

plot(y_2d_shift(:, 1, 1), C_trough_sp_2d, '--b');
h2=plot(y_2d_shift(1:3:end, 1, 1), C_trough_sp_2d(1:3:end), 'b^');

plot(y_3d_shift(:, 1, 1), C_peak_sp_3d, '-.g');
hold on;
h3=plot(y_3d_shift(1:3:end, 1, 1), C_peak_sp_3d(1:3:end), 'gp');

plot(y_3d_shift(:, 1, 1), C_trough_sp_3d, '--k');
h4=plot(y_3d_shift(1:3:end, 1, 1), C_trough_sp_3d(1:3:end), 'ks');

title('Spanwise', 'interpreter', 'latex');
hleg=legend([h1 h2 h3 h4],{'peak, $Ap/\delta = 2$','trough, $Ap/\delta = 2$', 'peak, $Ap/\delta = 3$','trough, $Ap/\delta = 3$'});
set(hleg, 'interpreter', 'latex');
xlim([0 5]);
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\rho_x(\Delta x, \Delta z = 0)$', 'interpreter', 'latex');
grid on;
grid minor;

%% 2D
x_4d_shift_2dimensional = repmat([fliplr(-1*x_4d_shift(:, :, 1)) x_4d_shift(:, 2:end, 1)], 2, 1);
x_4d_shift_2dimensional = x_4d_shift_2dimensional(1:199, :);

y_4d_shift_2dimensional = repmat(vertcat(flipud(y_4d_shift(2:end, :, 1)), -1*y_4d_shift(:, :, 1)), 1, 2);
y_4d_shift_2dimensional = y_4d_shift_2dimensional(:, 1:65);

figure()
set(gcf,'units','centimeters','position',[0 0 2.5*8 8]);

subplot(1,2,1)
contourf(x_4d_shift_2dimensional, y_4d_shift_2dimensional, C_mean_4d(:, :, 21)./C_mean_4d(100, 33, 21), 10);

colormap jet;
title('Local Maxima', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');
grid on;
grid minor;
caxis([-1 1]);

hb = colorbar('location','eastoutside');
w = hb.Position;
w(1) = 0.468;
w(2) = 0.14;
w(3) = 0.03;
w(4) = 0.79;
hb.Position = w;


subplot(1,2,2)
contourf(x_4d_shift_2dimensional, y_4d_shift_2dimensional, C_mean_4d(:, :, 45)./C_mean_4d(100, 33, 45), 10);

colormap jet;
colorbar;
title('Local Minima', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');
grid on;
grid minor;
caxis([-1 1]);

hb = colorbar('location','eastoutside');
w = hb.Position;
w(1) = 0.908;
w(2) = 0.14;
w(3) = 0.03;
w(4) = 0.79;
hb.Position = w;

% figure()
% set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 21)./C_mean_4d(100, 100, 21),'-.r');
% hold on;
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 76)./C_mean_4d(100, 100, 76), '-.r');
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 45)./C_mean_4d(100, 100, 45), '--b');
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 103)./C_mean_4d(100, 100, 103), '--b');
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 1)./C_mean_4d(100, 100, 1));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 10)./C_mean_4d(100, 100, 10));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 20)./C_mean_4d(100, 100, 20));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 30)./C_mean_4d(100, 100, 30));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 40)./C_mean_4d(100, 100, 40));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 50)./C_mean_4d(100, 100, 50));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 60)./C_mean_4d(100, 100, 60));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 70)./C_mean_4d(100, 100, 70));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 80)./C_mean_4d(100, 100, 80));
% plot(x_4d_shift(1, :, 1), C_mean_4d(100, 100:end, 90)./C_mean_4d(100, 100, 90));
% 
% plot(x_4d_shift(1, :, 1), C_peak_st_4d, '-.r');
% h1=plot(x_4d_shift(1, 1:3:end, 1), C_peak_st_4d(1:3:end), 'rd');
% plot(x_4d_shift(1, :, 1), C_trough_st_4d, '--b');
% h2=plot(x_4d_shift(1, 1:3:end, 1), C_trough_st_4d(1:3:end), 'b^');
% plot(x_4d_shift(1, :, 1), C_mean_complete_st, '-g');
% h3=plot(x_4d_shift(1, 1:3:end, 1), C_mean_complete_st(1:3:end), 'gx');
% 
% 
% xlim([0 12]);
% xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
% ylabel('$\rho_x(\Delta x, \Delta z = 0)$', 'interpreter', 'latex');
% grid on;
% grid minor;

%%
figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 21), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Maxima 1', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 21*2), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Maxima 2', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 21*3), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Maxima 3', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 76*1), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Maxima 4', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 76*2), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Maxima 5', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 76*3), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Maxima 6', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

%%%%%%%

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 45), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Minima 1', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 45*2), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Minima 2', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 45*3), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Minima 3', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 103*1), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Minima 4', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 103*2), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Minima 5', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');

figure();
set(gcf,'units','centimeters','position',[0 0 1.2*8 8]);
surf(x_4d_shift(:, :, 1), y_4d_shift(:, :, 1), phase_4d(:, :, 103*3), 'EdgeColor', 'none')
view(2);
colormap jet;
xlim([0 4]);
title('Local Minima 6', 'interpreter', 'latex');
xlabel('$\Delta x/\delta$', 'interpreter', 'latex');
ylabel('$\Delta y/\delta$', 'interpreter', 'latex');


