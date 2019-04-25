function [tbf] = tracbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,Xt,Yt,nt)
%[epsbf] = tracbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,Xt,Yt,nt)
%   INPUTS: lambda_Cx,lambda_Cy, [1,m] vectors containing lambda values for basis functions
%           L_Cx,L_Cy, scalars: basis function scaling
%           x,y, [n_*,1] vectors: test points
%   OUTPUT: d0bf [3*n_*,m] matrix: containing basis functions where each
%   [3,m] submatrix contains the basis fucntions for the strain function 
%   [eps_xx, eps_xy, eps_yy]^T


[dbf] = dbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,Xt(:),Yt(:));

dxdxC = dbf(:,:,1);
dydyC = dbf(:,:,2);
dxdyC = dbf(:,:,3);

n1 = nt(1,:)';
n2 = nt(2,:)';

[n,m,~] = size(dbf);
tbf = NaN(2*n,m);

tbf(1:2:2*n-1,:) = n1.*dydyC - n2.*dxdyC;
tbf(2:2:2*n,:) = -n1.*dxdyC + n2.*dxdxC;


end

