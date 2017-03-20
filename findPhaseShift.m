%%
clear all;
close all;

% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/piyushData2.mat');
load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed/F2_F14_processed.mat')
% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F2/Processed/F2_F06_processed.mat')
%%
N   = size(phase_Final, 2);
% trim = [5; 28; 6; 53];
trim = [4; 32; 3; 69];
% trim = [5; 27; 6; 53];


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
error = zeros(1, N);
start = 1;
for i = 1:N
    a    = phase_Final{start} - phase_Final{i};
    b    = a.^2;
    c    = sum(sum(b))./numel(b);
    rmse = sqrt(c);
    
    error(i) = rmse;
end

%%
C   = zeros(1, N);
x_c = 1:length(C);
for i = 1:N
    C(i) = sum(sum(phase_Final{start+60}.*phase_Final{i}, 1), 2);

end
[pks,locs] = findpeaks(C,'MinPeakDistance',40);
phaseShift = diff(locs);
mean_phaseShift = mean(phaseShift);
