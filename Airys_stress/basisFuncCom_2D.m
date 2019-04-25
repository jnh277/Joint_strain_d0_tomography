function [Idfuncs,lambdas,SLambda,dfuncs] = basisFuncCom_2D(MM,m,covFunc,lx,ly,sig_f,nhat,entry,exit,X,Y,nsegs)
% [Idfuncs,lambdas,SLambda,dfuncs] = basisFuncCom_2D(MM,m,covFunc,lx,ly,sig_f,nhat,entry,exit,X,Y,nsegs)
%   m: number of basis functions in each direction
%   MM: placing of basis functions to make a circular placement rather than
%   lx,ly: length scales
%   sig_f: prior uncertainty
%   nhat: ray directions, each column corresponds to a ray
%   entry: entry intersects for a ray, each column corresponds to a ray
%   exit: exit intersects for a ray, each column corresponds to a ray
%   X,Y: position of test points
%   nsegs: number of segments each ray intercepts

[~,n] = size(entry);        % number of ray measurements
[~,mm_adj] = size(MM);    

if ~exist('nsegs','var') || isempty(nsegs)
    nsegs = ones(1,n);
end


%% determine the problem domain and place the basis function lambdas
dlambda_x = 3.5/lx/m;           % 3.5 sigma coverage
dlambda_y = 3.5/ly/m;

lambda_x = MM(1,:)*dlambda_x;   % lambda spacings
lambda_y = MM(2,:)*dlambda_y;

Lx = pi/2/dlambda_x;            % Domain scales
Ly = pi/2/dlambda_y;


n1 = nhat(1,:)'; n2 = nhat(2,:)';

n1lx = bsxfun(@times,n1,lambda_x);
n2ly = bsxfun(@times,n2,lambda_y);

Q1 = (n1lx-n2ly);
Q2 = (n1lx+n2ly);

I1 = ceil(find(abs(Q1(:)) <= eps^(1/2))/length(n1));
I2 = ceil(find(abs(Q2(:)) <= eps^(1/2))/length(n1));

I = unique([I1;I2]);

% attempt to avoid divide by zero (or very small number) by moving the
% basis functions
while (~isempty(I))
    lambda_x(I) = lambda_x(I) + eps^(1/2)*max(lambda_x)*(1-rand(1,length(I)));
    lambda_y(I) = lambda_y(I) + eps^(1/2)*max(lambda_y)*(1-rand(1,length(I)));
    
    n1lx = bsxfun(@times,n1,lambda_x);
    n2ly = bsxfun(@times,n2,lambda_y);
    Q1 = (n1lx-n2ly);
    Q2 = (n1lx+n2ly);
    
    I1 = ceil(find(Q1(:) == 0)/length(n1));
    I2 = ceil(find(Q2(:) == 0)/length(n1));
    I = unique([I1;I2]);
end

% put the lambdas together for output
lambdas = [lambda_x;lambda_y];
% compute the spectral intensities of the basis functions
if strcmp('SE',covFunc)
    SLambda = sig_f^2*(2*pi)*lx*ly*exp(-0.5*(lambda_x.^2*lx^2+lambda_y.^2*ly^2));
elseif strcmp('M5_2',covFunc)
    SLambda = sig_f^2*4*pi*gamma(7/2)/gamma(5/2)*5^(2.5)*lx*ly./(5+lambda_x.^2*lx^2+lambda_y.^2*ly^2).^(3.5);
elseif strcmp('M3_2',covFunc)
    SLambda = sig_f^2*4*pi*gamma(5/2)/gamma(3/2)*3^(1.5)*lx*ly./(3+lambda_x.^2*lx^2+lambda_y.^2*ly^2).^(2.5);
elseif strcmp('M1_2',covFunc)
    SLambda = sig_f^2*4*pi*gamma(3/2)/gamma(1/2)*lx*ly./(1+lambda_x.^2*lx^2+lambda_y.^2*ly^2).^(1.5);
else
    error('Invalid covariance function')
end
    
    

% preallocate to zeros so multiseg stuff works  
Idfuncs = zeros(n,mm_adj,3);


%% Do the ray integrals of the function derivatives

for ss = 1:max(nsegs)
    segInd = find(nsegs >= ss);        % index of measurements that have at least this many segments
    
    x0 = entry(ss*2-1,segInd)'; y0 = entry(ss*2,segInd)';
    alpha_x = bsxfun(@times,Lx+x0,lambda_x);
    alpha_y = bsxfun(@times,Ly+y0,lambda_y);

    Gamma1s = sin(alpha_x-alpha_y);
    Gamma2s = sin(alpha_x+alpha_y);

    % exit points s = L
    xf = exit(ss*2-1,segInd)'; yf = exit(ss*2,segInd)';

    alpha_x = bsxfun(@times,Lx+xf,lambda_x);
    alpha_y = bsxfun(@times,Ly+yf,lambda_y);

    Gamma1f = sin(alpha_x-alpha_y);
    Gamma2f = sin(alpha_x+alpha_y);

    % merging
    Gamma1 = (Gamma1f-Gamma1s)./Q1(segInd,:);
    Gamma2 = (Gamma2f-Gamma2s)./Q2(segInd,:);
    
    % integral of d/dxdx
    Idfuncs(segInd,:,1) = Idfuncs(segInd,:,1) - (lambda_x.^2/sqrt(Lx*Ly)/2).*(Gamma1-Gamma2);
    
    % integral of d/dydy
    Idfuncs(segInd,:,2) = Idfuncs(segInd,:,2) -(lambda_y.^2/sqrt(Lx*Ly)/2).*(Gamma1-Gamma2);
    
    % integral of d/dxdy
    Idfuncs(segInd,:,3) = Idfuncs(segInd,:,3) +(lambda_x.*lambda_y/sqrt(Lx*Ly)/2).*(Gamma1+Gamma2);

end
%% Compute the basis function derivates
if nargout > 3
    n_test = length(X(:));
    dfuncs = NaN(n_test,mm_adj,3);
    
    Bx = bsxfun(@times,X(:)+Lx,lambda_x);       % here X,Y are for the test points
    By = bsxfun(@times,Y(:)+Ly,lambda_y);

    sx = sin(Bx);       % avoid calling the sin function lots
    sy = sin(By);
    cx = cos(Bx);
    cy = cos(By);

    % d/dxdx
    dfuncs(:,:,1) =  -(lambda_x.^2/sqrt(Lx*Ly)) .* sx .* sy;
    
    % d/dydy
    dfuncs(:,:,2) =  -(lambda_y.^2/sqrt(Lx*Ly)) .* sx .* sy;
    
    % d/dxdy
    dfuncs(:,:,3) =  ((lambda_x.*lambda_y)/sqrt(Lx*Ly)) .* cx .* cy;
    
    
end
    
end