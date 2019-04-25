function [LogL,grad] = hpOptim_airys_com(theta,num_basis,C,nu,entry,exit,y,nsegs,Pb,nb)
% [LogL,grad] = hpOptim_maxwell(l,num_basis,A,B,C,nu,entry,exit,y,sig_m,mode)

[~,n] = size(entry);
if length(theta) ~= 5
    error('Theta must contain sig_f, lx, ly, sig_t,sig_m')
end

sig_f = theta(1);
C.sig_f = sig_f;

lx = theta(2);
ly = theta(3);

C.lx = lx;
C.ly = ly;

sig_t = theta(4);

sig_m = theta(5);


if exist('nsegs','var') && ~isempty(nsegs)
    [Phi, SLambda,lambdas,Phi_bound] = airys_approx(num_basis,C,nu,entry,exit,Pb(:,1),Pb(:,2),nsegs);
else
    [Phi, SLambda,lambdas,Phi_bound] = airys_approx(num_basis,C,nu,entry,exit,Pb(:,1),Pb(:,2));
end

% ------------Include tractions--------------------------------
HookesS = 1/(1-nu^2)*[1, 0, nu;
                    0, 1-nu, 0;
                    nu, 0, 1]; 

HCell = repmat({HookesS}, 1, length(Pb));
HH = blkdiag(HCell{:}); 
Phi_bstress = HH*Phi_bound;
for i = 1:length(Pb)
    nbCell{i} = [nb(1,i), nb(2,i), 0;0,nb(1,i), nb(2,i)];
end
NbNb = blkdiag(nbCell{:});
Phi_trac =NbNb*Phi_bstress;

yTrac = zeros(length(Pb)*2,1);
Phi = [Phi;Phi_trac];
y = [y;yTrac];
nt = length(yTrac);
% ------------------------------------------------------

isnb = [ones(n,1)/sig_m^2;ones(nt,1)/sig_t^2];

[~,m] = size(Phi);          % total number of basies 
nm = m;            % number of basies per basis function (C)

% solving using QR instead for numerical fuckign reasons
Gamma = [[(ones(n,1)/sig_m);ones(nt,1)/sig_t].*Phi;diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);

logdetZ = 2*sum(log(abs(diag(CZ))));
logQ = sum(log(SLambda)) + n*log(sig_m^2) + logdetZ;

optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
v = (linsolve(CZ,linsolve(CZ,Phi'*(isnb.*y),optsT),opts));

yinvQy = y'*(isnb.*y) - (y.*isnb)'*Phi*v;

LogL = logQ/2 + yinvQy/2 + n/2*log(2*pi) + nt/2*log(2*pi);



if nargout == 2      
        % GRADIENTS OF C
        inds = 1:nm;
%         dlogQdsigf = 2/C.sig_f - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag(2/C.sig_f./SLambda(inds)))));
        dlogQdsigf = 2/C.sig_f - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag(2/C.sig_f./SLambda(inds)),optsT),opts));         
        dyinvQydsigf = -v(inds)'*diag(2/C.sig_f./SLambda(inds))*v(inds);
        grad_sigf = dlogQdsigf/2 + dyinvQydsigf/2;
        
%         dlogQdlx = sum(1/C.lx - C.lx*lambdas(1,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/C.lx - C.lx*lambdas(1,inds))./SLambda(inds)))));
        dlogQdlx = sum(1/C.lx - C.lx*lambdas(1,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/C.lx - C.lx*lambdas(1,inds))./SLambda(inds)),optsT),opts)); 
        dyinvQydlx = -v(inds)'*diag((1/C.lx - C.lx*lambdas(1,inds))./SLambda(inds))*v(inds);
        grad_lx = dlogQdlx/2 + dyinvQydlx/2;
        
%         dlogQdly = sum(1/C.ly - C.ly*lambdas(2,inds)) - trace(CZ(inds,inds)\(CZ(inds,inds)'\(diag((1/C.ly - C.ly*lambdas(2,inds))./SLambda(inds)))));
        dlogQdly = sum(1/C.ly - C.ly*lambdas(2,inds)) - trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),diag((1/C.ly - C.ly*lambdas(2,inds))./SLambda(inds)),optsT),opts)); 
        dyinvQydly = -v(inds)'*diag((1/C.ly - C.ly*lambdas(2,inds))./SLambda(inds))*v(inds);
        grad_ly = dlogQdly/2 + dyinvQydly/2;
        
        % GRADIENT with respect to TRACTION std
        vt = (linsolve(CZ,linsolve(CZ,Phi_trac'*(yTrac*2/sig_t^2),optsT),opts));
        dlogQdsigt = 2*nt/sig_t + trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),Phi_trac'*Phi_trac*2/sig_t^3,optsT),opts));
        dyinvQydsigt = - yTrac'*yTrac*2/sig_t^3 + 2*yTrac'*(2/sig_t^3)*Phi_trac*vt...
                        - vt'*(Phi_trac'*Phi_trac*2/sig_t^3)*vt;
        grad_sigt = dlogQdsigt/2 + dyinvQydsigt/2;
        
        % GRADIENT with respect to SIG_M
        vm = (linsolve(CZ,linsolve(CZ,Phi(1:n,:)'*(y(1:n)*2/sig_m^2),optsT),opts));
        dlogQdsigm = 2*n/sig_m + trace(linsolve(CZ(inds,inds),linsolve(CZ(inds,inds),Phi(1:n,:)'*Phi(1:n,:)*2/sig_m^3,optsT),opts));
        dyinvQydsigm = - y(1:n)'*y(1:n)*2/sig_m^3 + 2*y(1:n)'*(2/sig_m^3)*Phi(1:n,:)*vm...
                        - vm'*(Phi(1:n,:)'*Phi(1:n,:)*2/sig_m^3)*vm;
        grad_sigm = dlogQdsigm/2 + dyinvQydsigm/2;
        
        grad = [grad_sigf;grad_lx;grad_ly;grad_sigt;grad_sigm];
        


end
    


end