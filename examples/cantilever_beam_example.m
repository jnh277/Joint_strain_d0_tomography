clear all
clc

addpath ../helper_functions
addpath ../Airys_stress
addpath ../bayevar_airys

%% set whether to use already found hyperparameters or whether to run optimisation
% note turning optimisation on will result in the simulation taking quite a
% long time to complete, i.e. more than 10 hours
run_optimisation = false;    

%% simple cantilever beam strain field

% define parameters of cantilever beam
t = 6e-3;               % 6mm thickness             
h = 10e-3;              % 10mm height
l = 20e-3;              % 20mm length
E = 200e9;              % 200GPa young's modulus
nu = 0.3;                % poissons ratio
P = 2e3;                % 2KN load          % Is this meant to be positive or negative
I = t*h^3/12;           % second moment of inertia?

Fxx = @(x,y) P/(E*I).*(l-x).*y;               % idealised E_xx strain field
Fyy = @(x,y) -nu*P/(E*I).*(l-x).*y;           % idealised E_xy strain field
Fxy = @(x,y) -(1+nu)*P/(2*E*I)*(h^2/4 - y.^2);  % idealised E_yy strain field

bound = [0 h/2;l h/2;l -h/2;0 -h/2;0 h/2];
FV.faces = [1 2;
            2 3;
            3 4;
            4 1];


%% boundary and mesh


xb = bound(:,1);
yb = bound(:,2);

nx_t = 40;
ny_t = nx_t;
[X,Y] = meshgrid(linspace(0,l,nx_t),linspace(-h/2,h/2,ny_t));


IN = inpolygon(X,Y,xb,yb);

epsxx = nan(size(X));
epsxy = epsxx;
epsyy = epsxx;

epsxx(IN) = Fxx(X(IN),Y(IN));
epsxy(IN) = Fxy(X(IN),Y(IN));
epsyy(IN) = Fyy(X(IN),Y(IN));

max_epsxx = max(abs(epsxx(:)));
max_epsxy = max(abs(epsxy(:)));
max_epsyy = max(abs(epsyy(:)));


%% define a d0 variation field
mux = 0;
muy = 7e-3;
Lx = 15e-3/2;
Ly = 6e-3;
A = 0.0168;
d0_const = 4.056;

d0_field = @(x,y) A.*exp(-0.5*(x - mux).^2/Lx^2 - 0.5*(y-muy).^2/Ly^2) + d0_const;
d0_f = d0_field(X,Y);


%% traction stuff
numT = 50;
xt = [linspace(0,l,numT) linspace(l,0,numT)];
yt = [h/2*ones(1,numT) -h/2*ones(1,numT)];
Pb = [xt',yt'];
ntrac = [repmat([0;1],1,numT) repmat([0;-1],1,numT)];
%% create some measurement geometry
n_pixels = 100;
n_angles = 30;


angles = linspace(0,pi*(n_angles-1)/n_angles,n_angles);

[ pixelPoints ] = as_createPixelGrid(n_pixels,32*1e-3);
[ lines ] =  createRayLines(angles, pixelPoints);
lines = lines + [l/2, l/2;
                0, 0];

[mStructs, numM, shitlist,ylistr] = find_intersects(lines, bound(:,1),bound(:,2), 35e-3,4);
entry = [mStructs.entry];
exit = [mStructs.exit];
nhat = [mStructs.nHat];
n_obs = numel(mStructs);
L = [mStructs.L];

%% simulate the measurements using line integral
rng(1203)
sig_m = 1e-4;       % measurement variance in strain
[ dm ] = make_measurements_d0varying( entry, exit,Fxx, Fxy, Fyy,d0_field)+ d0_const*sig_m*randn(size(entry,2),1);
y_const = (dm - d0_const)/d0_const; % compute strain measurements assuming a constant d0 for the linear GP



%% set some parameters
sig_t = 1e-6;           % traction measurement variance
num_basis2 = 10;        % the number of basis functions in each x and y for the strain field
num_basis1 = 15;        % the number of basis functions in each x and y for the d0 field



if run_optimisation % warning this will take a while
    theta0 = [5;0.02;0.05;150;0.01;0.3];
    func = @(theta) var_inference_exp(theta,sig_m,dm,entry,exit,nhat,L,nu,d0_const,num_basis1,num_basis2,X(IN),Y(IN),xt(:),yt(:),ntrac,sig_t,true);
    optimSearch_options = optimset('display','iter','TolFun',1e-12,'TolX',1e-12,'MaxFunEvals',100,'PlotFcns',@optimplotfval);
    [thetaopt,Fval] = fminsearch(func,log(theta0),optimSearch_options);
    thetaopt = exp(thetaopt)

else    % use already found values (found values may differ depending on starting conditions due to local minima)
    theta_C = [182.4, 0.007, 0.2918]';
    theta_d0 = [6.6359, 0.0223, 0.0044]';
    thetaopt = [theta_d0;theta_C];
end


%% reconstruct using variational inference
disp('Maximising the free energy...')
[negF,wopt,Ropt,epsbf,d0bf] = ...
    var_inference(thetaopt,sig_m,dm,entry,exit,nhat,L,nu,d0_const,num_basis1,num_basis2,X(IN),Y(IN),xt(:),yt(:),ntrac,sig_t,true);
m2 = size(epsbf,2);
m1 = size(d0bf,2);
f = epsbf*wopt(m1+1:m1+m2);
epsxx_vi = nan(size(X));
epsxy_vi = nan(size(X));
epsyy_vi = nan(size(X));
epsxx_vi(:) = f(1:3:end);
epsxy_vi(:) = f(2:3:end);
epsyy_vi(:) = f(3:3:end);

d0_vi = nan(size(X));
d0_vi(:) = d0bf*wopt(1:m1);
%% reconstruct assuming constant d0 and using a linear GP method
C.sig_f = 150;  % pre optimised hyperparameters 
C.lx = 0.01;
C.ly = 0.3;
C.covFunc = 'SE';
P_all = [X(IN),Y(IN)];
P_new = [P_all;[xt(:), yt(:)]];     % include the boundary points for traction
[Phi, SLambda,lambdas,Phi_T] = airys_approx(num_basis2,C,nu,entry,exit,P_new(:,1),P_new(:,2),[]);
Phi_bound = Phi_T(length(P_all)*3+1:end,:);
Phi_pred = Phi_T(1:length(P_all)*3,:);
%% calculate the traction Phi
HookesS = 1/(1-nu^2)*[1, 0, nu;
                    0, 1-nu, 0;
                    nu, 0, 1]; 
HCell = repmat({HookesS}, 1, length(xt));
HH = blkdiag(HCell{:}); 
Phi_bstress = HH*Phi_bound;
for i = 1:length(xt)
    nbCell{i} = [ntrac(1,i), ntrac(2,i), 0;0,ntrac(1,i), ntrac(2,i)];
end
NbNb = blkdiag(nbCell{:});
Phi_trac =NbNb*Phi_bstress;

yTrac = zeros(length(xt)*2,1);
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

epsxx_pred = nan(size(X));
epsxy_pred = epsxx_pred;
epsyy_pred = epsxx_pred;

epsxx_pred(IN) = f_approx(1:3:end);
epsxy_pred(IN) = f_approx(2:3:end);
epsyy_pred(IN) = f_approx(3:3:end);

%%
errorGP = mean(abs([epsxx(IN);epsxy(IN);epsyy(IN)]-[epsxx_pred(IN);epsxy_pred(IN);epsyy_pred(IN)]))/max(abs([epsxx(IN);epsxy(IN);epsyy(IN)]))
errorVI = mean(abs([epsxx(IN);epsxy(IN);epsyy(IN)]-[epsxx_vi(IN);epsxy_vi(IN);epsyy_vi(IN)]))/max(abs([epsxx(IN);epsxy(IN);epsyy(IN)]))


%% plot theoretical strain field and both reconstructions
tmp = load('RdYlBu_up.mat');
cmap = tmp.RdYlBu;

clims_pred = 'auto';
clims_var = 'auto';
FigHandle = figure(4);
mm = 1e3;
micro = 1e6;

max_xx = max(abs([epsxx(:)]));
max_yy = max(abs([epsyy(:)]));
max_xy = max([epsxy(:)]);
min_xy = min([epsxy(:)]);

figure(4)
% epsxx
subplot(3,4,1); pcolor(X,Y,epsxx*micro); view(0,90); colorbar; colormap(cmap);  caxis([-max_xx max_xx]*micro); 
ylabel('FEA','FontSize',20)
shading interp; axis image; title('$\epsilon_{xx}$','Interpreter','latex','FontSize',30);
subplot(3,4,5); pcolor(X,Y,epsxx_pred*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_xx max_xx]*micro); 
shading interp; axis image;
ylabel('const d0 GP','FontSize',20)
subplot(3,4,9); pcolor(X,Y,epsxx_vi*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_xx max_xx]*micro); 
shading interp; axis image;
ylabel('VI GP','FontSize',20)

% epsxy
subplot(3,4,2); pcolor(X*mm,Y*mm,epsxy*micro); view(0,90); colorbar; colormap(cmap); caxis([min_xy max_xy]*micro); 
shading interp; axis image; title('$\epsilon_{xy}$','Interpreter','latex','FontSize',30);
subplot(3,4,6); pcolor(X*mm,Y*mm,epsxy_pred*micro); view(0,90); colorbar; colormap(cmap); caxis([min_xy max_xy]*micro); 
shading interp; axis image;
subplot(3,4,10); pcolor(X*mm,Y*mm,epsxy_vi*micro); view(0,90); colorbar; colormap(cmap); caxis([min_xy max_xy]*micro); 
shading interp; axis image;

% epsyy
subplot(3,4,3); pcolor(X*mm,Y*mm,epsyy*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_yy max_yy]*micro); 
shading interp; axis image; title('$\epsilon_{yy}$','Interpreter','latex','FontSize',30);
subplot(3,4,7); pcolor(X*mm,Y*mm,epsyy_pred*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_yy max_yy]*micro); 
shading interp; axis image;
subplot(3,4,11); pcolor(X*mm,Y*mm,epsyy_vi*micro); view(0,90); colorbar; colormap(cmap); caxis([-max_yy max_yy]*micro); 
shading interp; axis image;


% d0
subplot(3,4,4)
pcolor(X,Y,d0_f);axis image;shading interp;colorbar;title('$d_0$','Interpreter','latex','FontSize',30);

subplot(3,4,12)
pcolor(X,Y,d0_vi);axis image;shading interp;colorbar;

%% saving results
% save('CB_results.mat','X','Y','epsxx','epsyy','epsxy','epsxx_pred','epsxy_pred','epsyy_pred','epsxx_vi','epsxy_vi','epsyy_vi','d0_f','d0_vi')










