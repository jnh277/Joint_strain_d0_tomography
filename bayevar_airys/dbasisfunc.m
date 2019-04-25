function [dbf] = dbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,X,Y)
%[dbf] = dbasisfunc(lambda_Cx,lambda_Cy,L_Cx,L_Cy,x,y)
%   INPUTS: lambda_Cx,lambda_Cy, [1,m] vectors containing lambda values for basis functions
%           L_Cx,L_Cy, scalars: basis function scaling
%           x,y, [n_*,1] vectors: test points
%   OUTPUT: d0bf [3*n_*,m] matrix: containing basis functions where each
%   [3,m] submatrix contains the basis fucntions for the strain function 
%   [eps_xx, eps_xy, eps_yy]^T

m = length(lambda_Cx);
n_test = length(X(:));
dbf = NaN(n_test,m,3);

Bx = bsxfun(@times,X(:)+L_Cx,lambda_Cx);       % here X,Y are for the test points
By = bsxfun(@times,Y(:)+L_Cy,lambda_Cy);

sx = sin(Bx);       % avoid calling the sin function lots
sy = sin(By);
cx = cos(Bx);
cy = cos(By);

% d/dxdx
dbf(:,:,1) =  -(lambda_Cx.^2/sqrt(L_Cx*L_Cy)) .* sx .* sy;

% d/dydy
dbf(:,:,2) =  -(lambda_Cy.^2/sqrt(L_Cx*L_Cy)) .* sx .* sy;

% d/dxdy
dbf(:,:,3) =  ((lambda_Cx.*lambda_Cy)/sqrt(L_Cx*L_Cy)) .* cx .* cy;
    
end