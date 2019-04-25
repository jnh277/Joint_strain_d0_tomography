function [d0bf] = d0basisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,x,y)
%[d0bf] = d0basisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,x,y)
%   INPUTS: lambda_d0x,lambda_d0y, [1,m] vectors containing lambda values for basis functions
%           L_d0x,L_d0y, scalars: basis function scaling
%           x,y, [n_*,1] vectors: test points
%   OUTPUT: d0bf [n_*,m] matrix: containing basis functions


d0bf = sin(lambda_d0x.*(x+L_d0x)).*sin(lambda_d0y.*(y+L_d0y))./sqrt(L_d0x*L_d0y);
end

