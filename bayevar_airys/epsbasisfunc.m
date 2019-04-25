function [epsbf] = epsbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,X,Y,nu)
%[epsbf] = epsbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,x,y)
%   INPUTS: lambda_Cx,lambda_Cy, [1,m] vectors containing lambda values for basis functions
%           L_Cx,L_Cy, scalars: basis function scaling
%           x,y, [n_*,1] vectors: test points
%   OUTPUT: d0bf [3*n_*,m] matrix: containing basis functions where each
%   [3,m] submatrix contains the basis fucntions for the strain function 
%   [eps_xx, eps_xy, eps_yy]^T


[dbf] = dbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,X,Y);

dxdxC = dbf(:,:,1);
dydyC = dbf(:,:,2);
dxdyC = dbf(:,:,3);

m = length(lambda_Cx);
np = length(X(:));
epsbf = NaN(np*3,m);
epsbf(1:3:end,:) = dydyC-nu*dxdxC; % phi_exx;
epsbf(3:3:end,:) = dxdxC-nu*dydyC; % phi_eyy;
epsbf(2:3:end,:) = -(1+nu)*dxdyC;  % phi_exy;
end

