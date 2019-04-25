function [negF,wopt,Ropt,epsbf,d0bf] = ...
    var_inference_ms_exp(theta,sig_m,dm,entry,exit,nhat,nsegs,L,nu,d0_const,num_basis1,num_basis2,X,Y,Xt,Yt,ntrac,sig_t,display)
% Finds the approximate posterior that maximises the free energy
% this version of the code takes the log of the hyperparameters and is used
% for optimisation so that unconstrained optimisation can be used while
% making sure the parameters stay positive
% this function works with non convex geometry

if ~exist('display')
    display = false;
end

theta = exp(theta);

%% inputs
sig_fd0 = theta(1);
ldx = theta(2);
ldy = theta(3);
sig_fC = theta(4);
lCx = theta(5);
lCy = theta(6);

n = length(dm);
%% set up basis functions

[lambda_d0x,lambda_d0y,SLambdad0,lambda_Cx,lambda_Cy,SLambdaC,L_d0x,L_d0y,L_Cx,L_Cy] = ...
    placeBases(nhat,ldx,ldy,sig_fd0,lCx,lCy,sig_fC,num_basis1,num_basis2);
d0bf = d0basisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,X(:),Y(:));
[id0bf] = intd0basisfunc_ms(lambda_d0x,lambda_d0y,L_d0x,L_d0y,entry,exit,nhat,nsegs);
[ipbf,ipbf_order2] = intprodbasisfunc_ms(lambda_d0x,lambda_d0y,L_d0x,L_d0y,lambda_Cx,lambda_Cy,L_Cx,L_Cy,entry,exit,nhat,nsegs);



m1 = length(lambda_d0x);
m2 = length(lambda_Cx);
%% starting guess (hot starting d0 params)
% initialise d0 field params to be a constant field

yc = d0_const+0*X(:);
sig_c = 1e-8;
Gamma = [(ones(length(yc),1)/sig_c).*d0bf;diag(1./sqrt(SLambdad0))];
R = triu(qr(Gamma));
CZ = R(1:m1,1:m1);
wd00 = CZ\(CZ'\(d0bf'*yc/sig_c^2));

wC0 = zeros(m2,1);  % initialise strain field params to be a zero field

%% artificial traction measurements and the model
[tbf] = tracbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,Xt,Yt,ntrac);
nt = size(tbf,1);
ytrac = zeros(nt,1);

y = [dm;ytrac];
%% prepare for variational inference
mu = [wd00;wC0];

SLambda = [SLambdad0,SLambdaC];
Sigma = diag([ones(n,1)*sig_m^2;ones(nt,1)*sig_t^2]);
Sp = diag(SLambda);
max_iters = 15;
what = nan(m1+m2,max_iters+1);
what(:,1) = mu; % starting guess
Fhat = nan(max_iters+1,1);
dmhat = measurement_model(ipbf,id0bf,wd00,wC0,nhat,L,nu);
trhat = tbf*wC0;
yhat = [dmhat;trhat];
Fhat(1) = -0.5*(y - yhat)'*(Sigma\(y - yhat))-0.5*((mu-what(:,1))'/Sp)*(mu-what(:,1));
alpha_save = 1;
for i = 1:max_iters          % until convergence of free energy
    [Jwd0,JwC] = measurement_model_J(ipbf,ipbf_order2,id0bf,what(1:m1,i),what(m1+1:m1+m2,i),nhat,L,nu);
    J = [Jwd0,JwC;zeros(nt,m1),tbf];
    
    e1 = (y(1:n)-yhat(1:n));
    N1N1L = nhat(1,:)'.^2./L(:);
    N2N2L = nhat(2,:)'.^2./L(:);
    N1N2L = 2*nhat(1,:)'.*nhat(2,:)'./L(:);
    d2yhatdw2 =  N1N1L.*(ipbf(:,:,2)+nu*ipbf(:,:,1))-N1N2L.*(1+nu).*ipbf(:,:,3)+N2N2L.*(ipbf(:,:,1)+nu*ipbf(:,:,2));
    aH12 = reshape(sum((Sigma(1:n,1:n)\e1).*d2yhatdw2),m1,m2);

    g = J'*(Sigma\(y-yhat))+ Sp\(what(:,i) - mu);
    H = J'/Sigma*J+ diag(1./diag(Sp)) - [zeros(m1,m1), aH12;aH12',zeros(m2,m2)]+diag(1./diag(Sp));

    % for numerical stability
    [U,S,V] = svd(H);
    r = sum(diag(S) > eps*max(diag(S)));
    Ubar = U(:,1:r);
    Sbar = S(1:r,1:r);
    Vbar = V(:,1:r);
    Hg = Vbar*(Sbar\(Ubar'*g));
    

    Fnew_save = -inf;
    validPointFound = false;
    forward = true;     % initialy hope we get to increase step
    forward_steps = 0;      % steps attempted after a valid point is found
    backward_steps = 0;
    alpha = alpha_save;
    for k = 0:52        % a line search taht starts at the previously used step length 
        wnew = what(:,i)*(1-alpha) + alpha*mu + alpha*Hg;

        dmhat= measurement_model(ipbf,id0bf,wnew(1:m1),wnew(m1+1:m1+m2),nhat,L,nu);
        trhat = tbf*wnew(m1+1:m1+m2);
        yhatnew = [dmhat;trhat];
        Fnew = -0.5*(y - yhatnew)'*(Sigma\(y - yhatnew))-0.5*((mu-wnew)'/Sp)*(mu-wnew);
        
        if Fnew > Fhat(i) || validPointFound
            validPointFound = true;
            if Fnew > Fnew_save
                 alpha_save = alpha;
                 Fhat(i+1) = Fnew;
                 what(:,i+1) = wnew;
                 Fnew_save = Fnew;
                 yhat = yhatnew;
            elseif forward_steps>1 || backward_steps > 0
                break;
            else
               forward = false; 
               alpha = alpha_save;
            end
            if forward
                if alpha == 1       % don't do a larger step
                    break
                end
                 alpha = 1.5*alpha;
                 forward_steps = forward_steps+1;
            else
                 alpha = alpha/2;
                 backward_steps = backward_steps+1;       
            end
        else
            forward = false;
            alpha = alpha/2;
        end
        
    end
    if k==52
        i = i-1;
    end
    if display 
        disp(['Iter: ' num2str(i) ' Fval: ' num2str(Fhat(i+1)) ' Step: ' num2str(alpha_save) ' Fdif: ' num2str(abs(Fhat(i+1) - Fhat(i))) ' fcount: ' num2str(k+1) ' step norm: ' num2str(norm((what(:,i+1)-what(:,i))))])
    end
        
    if abs(Fhat(i+1) - Fhat(i)) < 1e-8*abs(Fhat(i+1)) || k == 52
        break
    end
end
Gamma = [diag(1./sqrt(SLambda));diag([ones(n,1)/sig_m;ones(nt,1)])*J];        % 
Ropt = triu(qr(Gamma));
wopt = what(:,i+1);

%%

[epsbf] = epsbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,X,Y,nu);

if display
figure(3)
subplot 211
plot(dm)
hold on
plot(yhat(1:n))
hold off
legend('measured dspacings','predicted dspacings')
title('predicted vs measured d spacings')
subplot 212
plot(yhat(n+1:n+nt))
title('predicted tractions')
end

%%
invR = pinv(Ropt);
logdetSigma = sum(log(diag(Sigma)));
logDetC = 2*sum(log(diag(abs(invR))));
logDetSigP = sum(log(SLambda));
sq1 = (mu - what(:,i+1))'*(Sp\(mu - what(:,i+1)));
sq2 = (y - yhat)'*(Sigma\(y - yhat));
negF = 0.5*(logdetSigma-logDetC+logDetSigP+sq1+sq2);

end

