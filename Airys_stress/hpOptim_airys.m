function [LogL,grad] = hpOptim_airys(theta,num_basis,C,nu,entry,exit,y,sig_m,nsegs)
% [LogL,grad] = hpOptim_maxwell(l,num_basis,A,B,C,nu,entry,exit,y,sig_m,mode)

n = length(y);
if length(theta) ~= 3
    error('Theta must contain sig_f, lx, ly')
end

sig_f = theta(1);
C.sig_f = sig_f;

lx = theta(2);
ly = theta(3);

C.lx = lx;
C.ly = ly;


if exist('nsegs','var')
    [Phi, SLambda,lambdas] = airys_approx(num_basis,C,nu,entry,exit,[],[],nsegs);
else
    [Phi, SLambda,lambdas] = airys_approx(num_basis,C,nu,entry,exit);
end

isnb = ones(n,1)/sig_m^2;

[n,m] = size(Phi);          % total number of basies 
nm = m;            % number of basies per basis function (C)

% solving using QR instead for numerical fuckign reasons
Gamma = [(ones(n,1)/sig_m).*Phi;diag(1./sqrt(SLambda))];
R = triu(qr(Gamma));
CZ = R(1:m,1:m);

logdetZ = 2*sum(log(abs(diag(CZ)))); % no multiply by 2 here as using qr not chol
logQ = sum(log(SLambda)) + n*log(sig_m^2) + logdetZ;

optsT.TRANSA = true; optsT.UT = true; opts.TRANSA = false; opts.UT = true;
v = (linsolve(CZ,linsolve(CZ,Phi'*(isnb.*y),optsT),opts));

yinvQy = y'*(isnb.*y) - (y.*isnb)'*Phi*v;

LogL = logQ/2 + yinvQy/2 + n/2*log(2*pi);



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

        
        grad = [grad_sigf;grad_lx;grad_ly];
        


end
    


end