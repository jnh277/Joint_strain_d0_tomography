function [id0bf] = intd0basisfunc_ms(lambda_d0x,lambda_d0y,L_d0x,L_d0y,entry,exit,nhat,nsegs)
%[id0bf] = intd0basisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,x,y)
%   INPUTS: lambda_d0x,lambda_d0y, [1,m] vectors containing lambda values for basis functions
%           L_d0x,L_d0y, scalars: basis function scaling
%           entry [2,n] vectors: where each column contains [x,y]^T
%           location of ray entering sample
%           exit [2,n] vectors: where each column contains [x,y]^T
%           location of ray exiting sample
%   OUTPUT: id0bf [n_*,m] matrix: containing line integral of basis
%   function

omega1 = zeros(length(nsegs),length(lambda_d0x));
omega2 = omega1;
for ss = 1:max(nsegs)
    segInd = find(nsegs >= ss);        % index of measurements that have at least this many segments
  
    % at entry
    x0 = entry(ss*2-1,segInd)'; y0 = entry(ss*2,segInd)'; 
    
    alpha_d0x = bsxfun(@times,x0+L_d0x,lambda_d0x);
    alpha_d0y = bsxfun(@times,y0+L_d0y,lambda_d0y);

    omega10 = sin(alpha_d0x - alpha_d0y);
    omega20 = sin(alpha_d0x + alpha_d0y);

    % at entry
    xf = exit(ss*2-1,segInd)'; yf = exit(ss*2,segInd)'; 
    
    alpha_d0x = bsxfun(@times,xf+L_d0x,lambda_d0x);
    alpha_d0y = bsxfun(@times,yf+L_d0y,lambda_d0y);

    omega11 = sin(alpha_d0x - alpha_d0y);
    omega21 = sin(alpha_d0x + alpha_d0y);
    
    omega1(segInd,:) = omega1(segInd,:) + omega11 - omega10;
    omega2(segInd,:) = omega2(segInd,:) + omega21 - omega20;

end
n1 = nhat(1,:)'; n2 = nhat(2,:)';

n1lx = bsxfun(@times,n1,lambda_d0x);
n2ly = bsxfun(@times,n2,lambda_d0y);

id0bf = (omega1./(2*(n1lx - n2ly)) ...
    - omega2./(2*(n1lx + n2ly)))/sqrt(L_d0x*L_d0y);



end

