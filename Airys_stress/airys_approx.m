function [Phi_yI, SLambda,lambdas,Phi_T] = airys_approx(num_basis,C,nu,entry,exit,X,Y,nsegs,L)
% [Phi_yI, SLambda,lambdas,Phi_T] = airys_approx(num_basis,C,nu,entry,exit,X,Y,nsegs,L)
% num_basis: nominal basis function resolution in each direction
% C: structures containing parameters for the underlying basis functions
% function
% X,Y: test points
% nu: poissons ratio
% entry: ray entry points (each point is a column/or each column contains the nsegs(i) points)
% exit: ray exit points (each point is a column)
% nsegs: number of segments each ray intercepted

m = num_basis;

if ~exist('nsegs','var')
    nsegs = [];
end

if isfield(C,'covFunc')
    covFunc= C.covFunc;
    if (~strcmp('SE',covFunc) && ~strcmp('M5_2',covFunc) && ~strcmp('M3_2',covFunc) && ~strcmp('M1_2',covFunc))
        error('Invalid covariance function')
    end
else
    covFunc = 'SE';
end

% generate basis function points
[mm1,mm2] = meshgrid(1:m);   % grid of basis functions  
% insideCircle = sqrt(mm1.^2 + mm2.^2)/m <=1+eps;    % points inside a circle
insideCircle = sqrt(mm1.^2 + mm2.^2)/m < 1;
mm1 = mm1(insideCircle);
mm2 = mm2(insideCircle);
MM = [mm1'; mm2'];   % all the basis functions, basis functions to change across columns
[~,mm_adj] = size(MM);

if ~exist('L','var')
    if isempty(nsegs) || all(nsegs == 1)
        L = sqrt(sum((entry-exit).^2));
%         nhat = (exit-entry)./L;
    else
        % check size of nsegs
        [r c] = size(nsegs);
        if r > 1
            nsegs = nsegs';
        end
        intmp = entry;
        extmp = exit;
        nan_inds = isnan(intmp(:));
        intmp(nan_inds) = 0;
        extmp(nan_inds) = 0;
        L= sqrt(sum((extmp-intmp).^2));
    end
else
    [r,c] = size(L);
    if r > c
        L = L';
    end
end

nhat = (exit(1:2,:)-entry(1:2,:))./hypot(exit(1,:)-entry(1,:),exit(2,:)-entry(2,:));



%% calculate Phi_yI

if nargout < 4

    [Idfuncs,lambdasC,SC] = basisFuncCom_2D(MM,m,covFunc,C.lx,C.ly,C.sig_f,nhat,entry,exit,[],[],nsegs);

    IdxdxC = Idfuncs(:,:,1);
    IdydyC = Idfuncs(:,:,2);
    IdxdyC = Idfuncs(:,:,3);
    


else

    [Idfuncs,lambdasC,SC,dfuncs] = basisFuncCom_2D(MM,m,covFunc,C.lx,C.ly,C.sig_f,nhat,entry,exit,X,Y,nsegs);

    IdxdxC = Idfuncs(:,:,1);
    IdydyC = Idfuncs(:,:,2);
    IdxdyC = Idfuncs(:,:,3);
    dxdxC = dfuncs(:,:,1);
    dydyC = dfuncs(:,:,2);
    dxdyC = dfuncs(:,:,3);
    
end

N1N1 = nhat(1,:)'.^2;
N2N2 = nhat(2,:)'.^2;
N1N2 = 2*nhat(1,:)'.*nhat(2,:)';


% relating to the 'C' basis functions for the rays
phi_yC = N1N1.*(IdydyC-nu*IdxdxC)+N2N2.*(IdxdxC-nu*IdydyC)+N1N2.*(-(1+nu)*IdxdyC); 
Phi_yI =(phi_yC)./L';

if nargout > 3
    np = length(X(:));
    % and combined for speed
    % interweave the Phis, so that we have [exx,exy,eyy] for each
    % point
    Phi_T = NaN(np*3,mm_adj);
    Phi_T(1:3:end,:) = dydyC-nu*dxdxC; % phi_exx;
    Phi_T(3:3:end,:) = dxdxC-nu*dydyC; % phi_eyy;
    Phi_T(2:3:end,:) = -(1+nu)*dxdyC;  % phi_exy;

end


% spectral values

SLambda = SC;
lambdas = lambdasC;
end





