%%
% phaseShift = 40;
% phaseShift = 55;
phaseShift = 111;

%%
phase_Final_new = phase_Final(68:end);

N = size(phase_Final_new, 2);
N  = size(phase_Final_new, 2) - mod(N, phaseShift);
NN = N/phaseShift;
%%
OPDrms = zeros(phaseShift, NN);
for i = 1:phaseShift
    for k = 1:NN
        a      = phase_Final_new{i+phaseShift*(k-1)};
        b      = a.^2;
        c      = sum(sum(b))./numel(b);
        OPDrms(i, k) = sqrt(c);
    end
end

OPDrms2 = mean(OPDrms, 2);
figure(); plot(OPDrms2);

%%
meanOPD = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);
for i = 1:phaseShift
    a = phase_Final_new(i:phaseShift:N);
    a = cat(3,a{:}); 
    
    meanOPD(:, :, i) = mean(a, 3);
    
end

varOPD = zeros(size(phase_Final_new{1}, 1), size(phase_Final_new{1}, 2), phaseShift);
for i = 1:phaseShift
    a = phase_Final_new(i:phaseShift:N);
    a = cat(3,a{:}); 
    
    varOPD(:, :, i) = var(a, 0,3);
    
end


%%

% v = VideoWriter('meanOPD.mp4', 'MPEG-4');
% v.FrameRate = 10;  % Default 30
% v.Quality   = 100;    % Default 75
% 
% open(v);
% 
% 
% for k = 1:phaseShift
%    surf(x_Final{k},y_Final{k},meanOPD(:, :, k),'edgecolor','none');
%    view(2);
%    colormap jet;
%    title(['mean(OPD(x, y)), Frame # ' num2str(k)]);
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% 
% close(v);
% 
% v = VideoWriter('varOPD.mp4', 'MPEG-4');
% v.FrameRate = 10;  % Default 30
% v.Quality   = 100;    % Default 75
% 
% open(v);
% 
% 
% for k = 1:phaseShift
%    surf(x_Final{k},y_Final{k},varOPD(:, :, k),'edgecolor','none');
%    view(2);
%    colormap jet;
%    title(['var(OPD(x, y)), Frame # ' num2str(k)]);
%    frame = getframe(gcf);
%    writeVideo(v,frame);
% end
% 
% close(v);

%%
varOPD_phi = squeeze(mean(mean(varOPD, 2), 1));
figure();
plot(varOPD_phi);

figure
errorbar(1:phaseShift,OPDrms2,varOPD_phi)
