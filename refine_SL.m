%%
clearvars;
close all;

% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/F3/Processed/F2_F14_processed.mat')

if ispc
    preamble = 'C:\Users\mkemnetz\Documents\MATLAB\Research\Shear Layer\shear-layer\F3';
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
fine = 100;
% nframes = 2000;
nframes = 45000;
% parfor i = 1:size(phase_Final, 2) 
parfor i = 1:nframes
    l = size(phase_Final{i}, 2);
    w = size(phase_Final{i}, 1);
    
    [Xq,Yq] = meshgrid(linspace(1,l, fine), linspace(1,w, fine));
    [X,Y]   = meshgrid(1:l, 1:w);
    
    phase_Final{i} = interp2(X, Y, phase_Final{i}, Xq, Yq, 'cubic');
    x_Final{i}     = interp2(X, Y, x_Final{i}    , Xq, Yq, 'spline');
    y_Final{i}     = interp2(X, Y, y_Final{i}    , Xq, Yq, 'spline');
end
phase_Final = phase_Final(1:nframes);
phase_m     = cat(3, phase_Final{:}); 

x_Final   = x_Final(1:nframes);
x_Final_m = cat(3, x_Final{:}); 

y_Final   = y_Final(1:nframes);
y_Final_m = cat(3, y_Final{:}); 
