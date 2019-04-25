clear all
clc
%% This script simulates measurements to be used with the C-shape example

%% mesh
nx = 100; ny = 100;
xmin = -10e-3;   xmax = 8e-3;
ymin = -10e-3;   ymax = 10e-3;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
% if ~any(y)==0 || isempty(y); y=[y(y<0) 0 y(y>0)]; ny=ny+1; end % for the slice plot
% y_zero_index = find(y==0);
[X,Y] = meshgrid(x,y); Xp=X; Yp=Y;
R = sqrt(X.^2+Y.^2); THETA = atan2(Y,X);
indices = R<=20e-3/2 & R>=7e-3/2 & (THETA>=pi/4 | THETA<=-pi/4); % points inside the C

theta = linspace(pi/4,2*pi-pi/4,100)';
[xtmp,ytmp] = pol2cart(theta,10e-3);
xb = [xtmp];
yb = [ytmp];

theta = linspace(2*pi-pi/4,pi/4,100)';
[xtmp,ytmp] = pol2cart(theta,3.5e-3);
xb = [xtmp;xb;xtmp(1)];
yb = [ytmp;yb;ytmp(1)];


%% artificial d0 field
mux = 5e-3;
muy = 7.5e-3;
Lx = 15e-3/2;
Ly = 6e-3;
% A = 0.01;
A = 0.0168;
d0_const = 4.056;

d0_field = @(x,y) A.*exp(-0.5*(x - mux).^2/Lx^2 - 0.5*(y-muy).^2/Ly^2) + d0_const;


d0_f = nan(size(X));

d0_f(indices) = d0_field(X(indices),Y(indices));

figure(5)
clf
pcolor(X,Y,d0_f)
hold on
plot(xb,yb)
hold off
title('d0 function')
shading interp

%% make rays and find intersects
n_pixels = 180;
n_angles = 60;


angles = linspace(0,pi*(n_angles-1)/n_angles,n_angles);


[ pixelPoints ] = as_createPixelGrid(n_pixels,30*1e-3);
[ lines ] =  createRayLines(angles, pixelPoints);


[mStructs, numM, shitlist,ylistr] = find_intersects(lines, xb,yb, 35e-3,4);
n_obs = length(mStructs);
nsegs = [mStructs.nsegs];         % allocate space for #segments per measurement

entry = NaN(max(nsegs)*2,n_obs);
exit = NaN(max(nsegs)*2,n_obs);
                                    % counter of total number of segments
entrytmp = [mStructs.entry];
entry(:) = entrytmp(:);
exittmp = [mStructs.exit];
exit(:) = exittmp(:);
nhat = [mStructs.nHat];
L = [mStructs.L]';

xlines = [entry(1,:);
    exit(1,:);
    nan(size(entry(1,:)))];

xlines2 = [entry(3,nsegs==2);
    exit(3,nsegs==2);
    nan(size(entry(3,nsegs==2)))];

ylines = [entry(2,:);
    exit(2,:);
    nan(size(entry(2,:)))];

ylines2 = [entry(4,nsegs==2);
    exit(4,nsegs==2);
    nan(size(entry(4,nsegs==2)))];

xlines = [xlines(:);xlines2(:)];
ylines = [ylines(:);ylines2(:)];



%%
load('C_shape_FEA.mat')

Exx = nan(size(X));
Exy = Exx;
Eyy = Exx;

Exx(indices) = Fxx(X(indices),Y(indices));
Exy(indices) = Fxy(X(indices),Y(indices));
Eyy(indices) = Fyy(X(indices),Y(indices));


figure(1)
clf

subplot 131
pcolor(X,Y,Exx)
shading interp
caxis([-max(abs(Exx(:))) max(abs(Exx(:)))])
axis equal
hold on
plot(xb,yb,'k','LineWidth',1.5)
hold off
title('Exx')

subplot 132
pcolor(X,Y,Exy)
shading interp
axis equal
hold on
plot(xb,yb,'k','LineWidth',1.5)
hold off
title('Exy')
caxis([-max(abs(Exy(:))) max(abs(Exy(:)))])

subplot 133
pcolor(X,Y,Eyy)
shading interp
axis equal
hold on
plot(xb,yb,'k','LineWidth',1.5)
hold off
title('Eyy')
caxis([-max(abs(Eyy(:))) max(abs(Eyy(:)))])


%% d0 varying measurements
sig_m = 1e-4;


d = zeros(n_obs,1);
for ss = 1:max(nsegs)
    segs = nsegs >= ss;
    Lseg = sqrt(sum((exit(ss*2-1:ss*2,segs)-entry(ss*2-1:ss*2,segs)).^2,1)).';
    d(segs) = d(segs) + make_measurements_d0varying( entry(ss*2-1:ss*2,segs), exit(ss*2-1:ss*2,segs),Fxx, Fxy, Fyy,d0_field).*Lseg;
end
d = d./L;

dm = d  + d0_const*sig_m*randn(size(entry,2),1);
y_const = (dm - d0_const)/d0_const; % measurement assuming constant d0

%% save
save('C_shape_sim_data3_60','entry','exit','nhat','nsegs','y_const','d','dm','L')


