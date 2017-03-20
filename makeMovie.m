% clear all;
% close all;
filelead = 'test';
% load('/Volumes/LaCie/MATLAB/Research/Shear Layer/piyushData2.mat')

loops = 500;
%%

% Apply the mask to the WFs
% F01 = (8:39, 7:40)  or trim = [8; 39; 7; 40];
% F01 = (11:36, 3:43) or trim = [11; 36; 3; 43];

% trim = [6; 27; 3; 57];
trim = [2; 21; 2; 21];

trim_a = trim(1);
trim_b = trim(2);
trim_c = trim(3);
trim_d = trim(4);
parfor i = 1:loops+5;
    phase_Final{i} = ...
        phase_Final{i}(trim_a:trim_b, trim_c:trim_d);
    x_Final{i}     = ...
        x_Final{i}(trim_a:trim_b, trim_c:trim_d);
    y_Final{i}     = ...
        y_Final{i}(trim_a:trim_b, trim_c:trim_d); 
end

%  Set the number of rows and columns for the interpolated wavefront.  The
%  main algorithm uses interp2 to refine the mesh.  
fine = 200;
parfor i = 1:loops+5;
    l = size(phase_Final{i}, 2);
    w = size(phase_Final{i}, 1);

    [Xq,Yq] = meshgrid(linspace(1,l, fine), linspace(1,w, fine));
    [X,Y]   = meshgrid(1:l, 1:w);

    phase_Final{i} = interp2(X, Y, phase_Final{i}, Xq, Yq, 'cubic');
    x_Final{i}     = interp2(X, Y, x_Final{i}    , Xq, Yq, 'spline');
    y_Final{i}     = interp2(X, Y, y_Final{i}    , Xq, Yq, 'spline');
end

%%


loops = 500;
Mov(loops) = struct('cdata',[],'colormap',[]);


for iframe=1:loops

surf(x_Final{iframe},y_Final{iframe},phase_Final{iframe},'edgecolor','none')



shading interp;
   colormap jet;
   box off;
   set(gcf, 'color', 'w');   
%    axis([0 2*pi 0 2*pi -400 400]);
    
   axis equal;
%    zlim([-1 1]);
   title(['WF(x, y), Frame # ' num2str(iframe)]);
   
   
   %You could use something like this to change the view
%    count = count+1;
%    view(-count,45)
%    
   %Or zoom in and out?
   %zoom(2*cos(i/25) - 0.2)
   
   %This is the stuff that makes the surface pretty, need to play with it
%    lightangle(-45,50)
%    h.FaceLighting = 'gouraud';
%    h.AmbientStrength = 0.50;
%    h.DiffuseStrength = 0.875;
%    h.SpecularStrength = 0.05;
%    h.SpecularExponent = 25;
%    h.BackFaceLighting = 'unlit';
%     view(0,75)
    shading interp
    lightangle(-45,30)
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.9;
    h.SpecularExponent = 25;
    h.BackFaceLighting = 'unlit';
   
   Mov(iframe)=getframe(gcf);

   
end

movie2avi(Mov,['Movie',filelead,'.avi'],'fps',10)
%%
%# preallocate
nFrames = loops;
f = getframe(gca);
[f,map] = rgb2ind(f.cdata, 256, 'nodither');
mov = repmat(f, [1 1 1 nFrames]);

%# create movie
for k=1:nFrames
    surf(x_Final{k},y_Final{k},phase_Final{k},'edgecolor','none')
    
    shading interp;
   colormap jet;
   box off;
   set(gcf, 'color', 'w');   
%    axis([0 2*pi 0 2*pi -400 400]);
    
   axis equal;
%    zlim([-1 1]);
   title(['WF(x, y), Frame # ' num2str(k)]);
   
   
   %You could use something like this to change the view
%    count = count+1;
%    view(-count,45)
%    
   %Or zoom in and out?
   %zoom(2*cos(i/25) - 0.2)
   
   %This is the stuff that makes the surface pretty, need to play with it
%    lightangle(-45,50)
%    h.FaceLighting = 'gouraud';
%    h.AmbientStrength = 0.50;
%    h.DiffuseStrength = 0.875;
%    h.SpecularStrength = 0.05;
%    h.SpecularExponent = 25;
%    h.BackFaceLighting = 'unlit';
%     view(0,75)
    shading interp
    lightangle(-45,30)
    h.FaceLighting = 'gouraud';
    h.AmbientStrength = 0.3;
    h.DiffuseStrength = 0.8;
    h.SpecularStrength = 0.9;
    h.SpecularExponent = 25;
    h.BackFaceLighting = 'unlit';
    
    
    
    f = getframe(gca);
    mov(:,:,1,k) = rgb2ind(f.cdata, map, 'nodither');
end
close(gcf)

%# create GIF and open
imwrite(mov, map, 'myPeaks4.gif', 'DelayTime',0, 'LoopCount',inf)
