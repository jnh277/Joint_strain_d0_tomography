function [LogL,grad] = hpOptim_airys_d0(theta,num_basis,C,nu,entry,exit,y,sig_m,nsegs,Pb,nb,sig_t,L)
% [LogL,grad] = hpOptim_maxwell(l,num_basis,A,B,C,nu,entry,exit,y,sig_m,mode)

n = length(y);
L = L(:);   % make sure L is a column vector
ndc = length(theta);
d0c = theta;

% sig_f = C.sig_f;
% lx = C.lx;
% ly = C.ly;

if exist('nsegs','var') && ~isempty(nsegs)
    [Phi, SLambda,lambdas,Phi_bound] = airys_approx(num_basis,C,nu,entry,exit,Pb(:,1),Pb(:,2),nsegs);
else
    [Phi, SLambda,lambdas,Phi_bound] = airys_approx(num_basis,C,nu,entry,exit,Pb(:,1),Pb(:,2));
end

%% d0 map
dPhidd0 = Phi;

d0s = polyval(d0c,L);
d0Map = diag(d0s);
Phi = d0Map*Phi;

%% tractions
HookesS = 1/(1-nu^2)*[1, 0, nu;
                    0, 1-nu, 0;
                    nu, 0, 1]; 

HCell = repmat({HookesS}, 1, length(Pb));
HH = blkdiag(HCell{:}); 
Phi_bstress = HH*Phi_bound;
clear nbCell;
for i = 1:length(Pb)
    nbCell{i} = [nb(1,i), nb(2,i), 0;0,nb(1,i), nb(2,i)];
end
NbNb = blkdiag(nbCell{:});
Phi_trac =NbNb*Phi_bstress;

yTrac = zeros(length(Pb)*2,1);
Phi = [Phi;Phi_trac];

%%
y = [y-d0s;yTrac];
%%
[nt,~] = size(yTrac);

isnb = [ones(n,1)/sig_m^2;ones(nt,1)/sig_t^2];

[~,m] = size(Phi);          % total number of basies 
nm = m;                     % number of basies per basis function (C)

% solving using QR instead for numerical fuckign reasons
Gamma = [[ones(n,1)/sig_m;ones(nt,1)/sig_t].*Phi;diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);

logdetZ = 2*sum(log(abs(diag(CZ))));
logQ = sum(log(SLambda)) + n*log(sig_m^2) + logdetZ;

optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
v = (linsolve(CZ,linsolve(CZ,Phi'*(isnb.*y),optsT),opts));

yinvQy = y'*(isnb.*y) - (y.*isnb)'*Phi*v;

LogL = logQ/2 + yinvQy/2 + n/2*log(2*pi);


if nargout == 2    
    grad = NaN(ndc,1);
    inds = 1:nm;
    
    % gradient with respect to d0 coefficients
    for i = 1:ndc
        dlogQdd0 = trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),2*(L.^(ndc-1).*dPhidd0)'*(isnb(1:n).*Phi(1:n,:)),optsT),opts));
        dyinvQydd0 = -2*isnb(1:n)'*y(1:n) + 2*isnb(1:n)'*Phi(1:n,:)*v...
            - 2*(isnb(1:n).*y(1:n))'*(L.^(ndc-1).*dPhidd0)*v+2*v'*(L.^(ndc-1).*dPhidd0)'*(isnb(1:n).*Phi(1:n,:))*v;
        grad_d0i = dlogQdd0/2+dyinvQydd0/2;
        
        grad(i) = grad_d0i;
    end


end
    


end