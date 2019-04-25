function [id0bf] = intd0basisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,entry,exit,nhat)
%[id0bf] = intd0basisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,x,y)
%   INPUTS: lambda_d0x,lambda_d0y, [1,m] vectors containing lambda values for basis functions
%           L_d0x,L_d0y, scalars: basis function scaling
%           entry [2,n] vectors: where each column contains [x,y]^T
%           location of ray entering sample
%           exit [2,n] vectors: where each column contains [x,y]^T
%           location of ray exiting sample
%   OUTPUT: id0bf [n_*,m] matrix: containing line integral of basis
%   function

% at entry
alpha_d0x = bsxfun(@times,entry(1,:)'+L_d0x,lambda_d0x);
alpha_d0y = bsxfun(@times,entry(2,:)'+L_d0y,lambda_d0y);

omega10 = sin(alpha_d0x - alpha_d0y);
omega20 = sin(alpha_d0x + alpha_d0y);

% at entry
alpha_d0x = bsxfun(@times,exit(1,:)'+L_d0x,lambda_d0x);
alpha_d0y = bsxfun(@times,exit(2,:)'+L_d0y,lambda_d0y);

omega11 = sin(alpha_d0x - alpha_d0y);
omega21 = sin(alpha_d0x + alpha_d0y);

n1 = nhat(1,:)'; n2 = nhat(2,:)';

n1lx = bsxfun(@times,n1,lambda_d0x);
n2ly = bsxfun(@times,n2,lambda_d0y);

id0bf = ((omega11 - omega10)./(2*(n1lx - n2ly)) ...
    - (omega21  - omega20)./(2*(n1lx + n2ly)))/sqrt(L_d0x*L_d0y);



end

