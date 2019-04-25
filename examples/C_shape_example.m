%% Recon C shape
clc
clear all

addpath ../helper_functions
addpath ../Airys_stress
addpath ../bayevar_airys

%% set whether to use already found hyperparameters or whether to run optimisation
% note turning optimisation on will result in the simulation taking quite a
% long time to complete, i.e. more than 10 hours
run_optimisation = false;    

%% mesh
nx = 100; ny = 100;
xmin = -10e-3;   xmax = 8e-3;
ymin = -10e-3;   ymax = 10e-3;
x = linspace(xmin,xmax,nx);
y = linspace(ymin,ymax,ny);
[X,Y] = meshgrid(x,y); Xp=X; Yp=Y;
R = sqrt(X.^2+Y.^2); THETA = atan2(Y,X);
indices = R<=20e-3/2 & R>=7e-3/2 & (THETA>=pi/4 | THETA<=-pi/4); % points inside the C


theta = linspace(pi/4,2*pi-pi/4,100)';
[xtmp,ytmp] = pol2cart(theta,10e-3);
xb = [xtmp];
yb = [ytmp];

inds = ~(((theta > deg2rad(70)) .* (theta < deg2rad(110))) + ((theta > deg2rad(250)) .* (theta < deg2rad(290))));
xt = xtmp(inds); % tractions
yt = ytmp(inds);
ntrac = [cos(theta(inds).');sin(theta(inds).')];


theta = linspace(2*pi-pi/4,pi/4,100)';
[xtmp,ytmp] = pol2cart(theta,3.5e-3);
xb = [xtmp;xb;xtmp(1)];
yb = [ytmp;yb;ytmp(1)];

inds = 1:3:length(xtmp);
xt = [xt;xtmp(inds)];
yt = [yt;ytmp(inds)];
ntrac = [ntrac, -[cos(theta(inds).');sin(theta(inds).')]];

xtmp = linspace(xb(end-1),xb(end),10)';
ytmp = linspace(yb(end-1),yb(end),10)';
theta = repmat(pi/4,10,1);

xt = [xt;xtmp];
yt = [yt;ytmp];
ntrac = [ntrac, [cos(theta.');sin(theta.')]];

xtmp = linspace(xb(100),xb(100+1),15)';
ytmp = linspace(yb(100),yb(100+1),15)';
theta = repmat(pi/4-pi/2,15,1);

xt = [xt;xtmp];
yt = [yt;ytmp];
ntrac = [ntrac, [cos(theta.');sin(theta.')]];

pc = [xt.';yt.'];

%% artificial d0 field
mux = 5e-3;
muy = 7.5e-3;
Lx = 15e-3/2;
Ly = 6e-3;
A = 0.01;
d0_const = 4.056;
d0_field = @(x,y) A.*exp(-0.5*(x - mux).^2/Lx^2 - 0.5*(y-muy).^2/Ly^2) + d0_const;


d0_f = nan(size(X));
d0_f(indices) = d0_field(X(indices),Y(indices));


%% load measurements 
load('C_shape_sim_data4_60.mat')
badind = isnan(d);
y_const(badind) = [];
d(badind) = [];
dm(badind) = [];
entry(:,badind) = [];
exit(:,badind) = [];
nhat(:,badind) = [];
nsegs(:,badind) = [];
L(badind) = [];
%%
%% setting some parameters

num_basis2 = 15;        % the number of basis functions in each x and y for the strain field
num_basis1 = 13;        % the number of basis functions in each x and y for the d0 field


nu = 0.3;
sig_m = 1e-4;       % measurement variance
sig_t = 5e-6;       % traction variance

%% attempt to optimise parameters
% 
if run_optimisation % warning this will take a while
    theta0 = [6.6359, 0.0223, 0.0044, 100, 0.0057, 0.0057]';
    func = @(theta) var_inference_ms_exp(theta,sig_m,dm,entry,exit,nhat,nsegs,L,nu,d0_const,num_basis1,num_basis2,X(indices),Y(indices),xt,yt,ntrac,sig_t,true);
    optimSearch_options = optimset('display','iter','TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',100,'PlotFcns',@optimplotfval);
    [thetaopt,Fval] = fminsearch(func,log(Theta),optimSearch_options);

    Theta = exp(thetaopt)
else    % use prefound optimisation results
    Theta = [7.0071;
        0.0185;
        0.0046;
       50.0328;
        0.0035;
        0.0046];
end

%% reconstruct with variational inference


[negF,wopt,Ropt,epsbf,d0bf] = ...
    var_inference_ms(Theta,sig_m,dm,entry,exit,nhat,nsegs,L,nu,d0_const,num_basis1,num_basis2,X(indices),Y(indices),xt,yt,ntrac,sig_t,true);
m2 = size(epsbf,2);
m1 = size(d0bf,2);
f = epsbf*wopt(m1+1:m1+m2);

epsxx_vi = nan(size(X));
epsxy_vi = epsxx_vi;
epsyy_vi = epsxx_vi;
d0_vi = epsxx_vi;

epsxx_vi(indices) = f(1:3:end);
epsxy_vi(indices) = f(2:3:end);
epsyy_vi(indices) =  f(3:3:end);
d0_vi(indices) = d0bf*wopt(1:m1);


%% reconstruct assuming a constant d0
C.sig_f = 50.0328; % using the prefound optimal hyper parameters
C.lx = 0.0035;
C.ly = 0.0046;
C.covFunc = 'SE';
P_all = [X(indices),Y(indices)];

P_new = [P_all;[xt, yt]];     % include the boundary points for traction
[Phi, SLambda,lambdas,Phi_T] = airys_approx(num_basis2,C,nu,entry,exit,P_new(:,1),P_new(:,2),nsegs);
Phi_bound = Phi_T(length(P_all)*3+1:end,:);
Phi_pred = Phi_T(1:length(P_all)*3,:);

%% calculate the traction Phi
HookesS = 1/(1-nu^2)*[1, 0, nu;
                    0, 1-nu, 0;
                    nu, 0, 1]; 
HCell = repmat({HookesS}, 1, length(pc));
HH = blkdiag(HCell{:}); 
Phi_bstress = HH*Phi_bound;
for i = 1:length(pc)
    nbCell{i} = [ntrac(1,i), ntrac(2,i), 0;0,ntrac(1,i), ntrac(2,i)];
end
NbNb = blkdiag(nbCell{:});
Phi_trac =NbNb*Phi_bstress;

yTrac = zeros(length(pc)*2,1);
nt = length(yTrac);
%% combine things
Phi = [Phi;Phi_trac];


%%

ys = [y_const/sig_m^2;yTrac/sig_t^2];
[~,m]= size(Phi);
n = length(y_const);
Gamma = [[(ones(n,1)/sig_m);ones(nt,1)/sig_t].*Phi;diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);
optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
f_approx = Phi_pred*(linsolve(CZ,linsolve(CZ,Phi'*ys,optsT),opts));

epsxx_GP = nan(size(X));
epsxy_GP = epsxx_GP;
epsyy_GP = epsxx_GP;


epsxx_GP(indices) = f_approx(1:3:end);
epsxy_GP(indices) = f_approx(2:3:end);
epsyy_GP(indices) = f_approx(3:3:end);


%% load FEA
load('C_shape_FEA.mat')

Exx = nan(size(X));
Exy = Exx;
Eyy = Exx;

Exx(indices) = Fxx(X(indices),Y(indices));
Exy(indices) = Fxy(X(indices),Y(indices));
Eyy(indices) = Fyy(X(indices),Y(indices));

%% plot fea and both reconstructions reconstruction
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;

clims_pred = 'auto';
clims_var = 'auto';
FigHandle = figure(4);
mm = 1e3;
micro = 1e6;

max_xx = max(abs([Exx(:)]));
max_yy = max(abs([Eyy(:)]));
max_xy = max([Exy(:)]);
min_xy = min([Exy(:)]);

figure(4)
% epsxx
subplot(3,4,1); pcolor(Xp,Yp,Exx*micro); view(0,90); colorbar; colormap(cmap);  caxis([-max_xx max_xx]*micro); 
ylabel('FEA','FontSize',20)
shading interp; axis image; title('$\epsilon_{xx}$','Interpreter','latex','FontSize',30);
subplot(3,4,5); pcolor(Xp,Yp,epsxx_GP*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_xx max_xx]*micro); 
shading interp; axis image;
ylabel('const d0 GP','FontSize',20)

subplot(3,4,9); pcolor(Xp,Yp,epsxx_vi*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_xx max_xx]*micro); 
shading interp; axis image;
ylabel('VI GP','FontSize',20)

% epsxy
subplot(3,4,2); pcolor(Xp*mm,Yp*mm,Exy*micro); view(0,90); colorbar; colormap(cmap); caxis([min_xy max_xy]*micro); 
shading interp; axis image; title('$\epsilon_{xy}$','Interpreter','latex','FontSize',30);
subplot(3,4,6); pcolor(Xp*mm,Yp*mm,epsxy_GP*micro); view(0,90); colorbar; colormap(cmap); caxis([min_xy max_xy]*micro); 
shading interp; axis image;
subplot(3,4,10); pcolor(Xp*mm,Yp*mm,epsxy_vi*micro); view(0,90); colorbar; colormap(cmap); caxis([min_xy max_xy]*micro); 
shading interp; axis image;

% epsyy
subplot(3,4,3); pcolor(Xp*mm,Yp*mm,Eyy*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_yy max_yy]*micro); 
shading interp; axis image; title('$\epsilon_{yy}$','Interpreter','latex','FontSize',30);
subplot(3,4,7); pcolor(Xp*mm,Yp*mm,epsyy_GP*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_yy max_yy]*micro); 
shading interp; axis image;
subplot(3,4,11); pcolor(Xp*mm,Yp*mm,epsyy_vi*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_yy max_yy]*micro); 
shading interp; axis image;

% d0
d0_f = nan(size(Xp));
d0_f(indices) = d0_field(Xp(indices),Yp(indices));
subplot(3,4,4)
pcolor(Xp,Yp,d0_f);axis image;shading interp;colorbar;title('$d_0$','Interpreter','latex','FontSize',30);

subplot(3,4,12)
pcolor(Xp,Yp,d0_vi);axis image;shading interp;colorbar;

%%
errorGP = mean(abs([Exx(indices);Exy(indices);Eyy(indices)]-[epsxx_GP(indices);epsxy_GP(indices);epsyy_GP(indices)]))/max(abs([Exx(indices);Exy(indices);Eyy(indices)]))

errorVI = mean(abs([Exx(indices);Exy(indices);Eyy(indices)]-[epsxx_vi(indices);epsxy_vi(indices);epsyy_vi(indices)]))/max(abs([Exx(indices);Exy(indices);Eyy(indices)]))


