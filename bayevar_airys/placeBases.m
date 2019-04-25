function [lambda_d0x,lambda_d0y,SLambdad0,lambda_Cx,lambda_Cy,SLambdaC,L_d0x,L_d0y,L_Cx,L_Cy] = ...
    placeBases(nhat,ldx,ldy,sig_fd0,lCx,lCy,sig_fC,num_basis1,num_basis2)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
m1 = num_basis1;

[mm1,mm2] = meshgrid(1:m1);   % grid of basis functions  
insideCircle = sqrt(mm1.^2 + mm2.^2)/m1 < 1;
mm1 = mm1(insideCircle);
mm2 = mm2(insideCircle);
MM = [mm1'; mm2'];   % all the basis functions, basis functions to change across columns
% [~,mm_adj1] = size(MM);

dlambda_x = 4.5/ldx/m1;           % 4 sigma coverage
dlambda_y = 4.5/ldy/m1;
lambda_d0x = MM(1,:)*dlambda_x;   % lambda spacings
lambda_d0y = MM(2,:)*dlambda_y;
L_d0x = pi/2/dlambda_x;            % Domain scales
L_d0y = pi/2/dlambda_y;


m2 = num_basis2;

[mm1,mm2] = meshgrid(1:m2);   % grid of basis functions  
insideCircle = sqrt(mm1.^2 + mm2.^2)/m2 < 1;
mm1 = mm1(insideCircle);
mm2 = mm2(insideCircle);
MM = [mm1'; mm2'];   % all the basis functions, basis functions to change across columns

dlambda_x = 3.5/lCx/m2;           % 3.5 sigma coverage
dlambda_y = 3.5/lCy/m2;
lambda_Cx = MM(1,:)*dlambda_x;   % lambda spacings
lambda_Cy = MM(2,:)*dlambda_y;
L_Cx = pi/2/dlambda_x;            % Domain scales
L_Cy = pi/2/dlambda_y;


%% This entire block is for protection against divide by zeros
n1 = nhat(1,:)';
n2 = nhat(2,:)';
n1ldx = bsxfun(@times,n1,lambda_d0x);
n2ldy = bsxfun(@times,n2,lambda_d0y);

m2 = length(lambda_Cx);
m1 = length(lambda_d0x);        % on purpose overwriting the old

Lamdx = repmat(lambda_d0x(:),1,m2); Lamdx = Lamdx(:)';
Lamdy = repmat(lambda_d0y(:),1,m2); Lamdy = Lamdy(:)'; % d0 changes in each sub block
LamCx = repmat(lambda_Cx(:)',m1,1); LamCx = LamCx(:)';
LamCy = repmat(lambda_Cy(:)',m1,1); LamCy = LamCy(:)';

n1Ldx = bsxfun(@times,n1,Lamdx);
n2Ldy = bsxfun(@times,n2,Lamdy);
n1LCx = bsxfun(@times,n1,LamCx);
n2LCy = bsxfun(@times,n2,LamCy);

Q1 = abs(n1Ldx-n1LCx+n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % - + +
Q2 = abs(n1Ldx-n1LCx-n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % - - -
Q3 = abs(n1Ldx-n1LCx+n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % - + -
Q4 = abs(n1Ldx+n1LCx-n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % + - - 
Q5 = abs(n1Ldx+n1LCx+n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % + + -
Q6 = abs(n1Ldx-n1LCx-n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % - - +
Q7 = abs(n1Ldx+n1LCx-n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % + - +
Q8 = abs(n1Ldx+n1LCx+n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % + + +

Q1d = abs(n1ldx-n2ldy) < eps^(2/3)*max(lambda_d0x);
Q2d = abs(n1ldx+n2ldy) < eps^(2/3)*max(lambda_d0x);


I = logical(sum(Q1d+Q2d));
I2 = logical(sum(Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8));

% shift = sqrt(eps);
maxVal = max([lambda_Cx, lambda_Cy, lambda_d0x, lambda_d0y]);
shift = 1e-6;
while any(I) || any(I2)
    lambda_d0x(I) = lambda_d0x(I) + shift*max(lambda_d0x)*(1-rand(1,sum(I)));
    lambda_d0y(I) = lambda_d0y(I) + shift*max(lambda_d0x)*(1-rand(1,sum(I)));
    
    I2Mat = nan(m1,m2);
    I2Mat(:) = I2;
    I2d = logical(sum(I2Mat,2))';
    I2C = logical(sum(I2Mat,1));
    lambda_d0x(I2d) = lambda_d0x(I2d) + shift*max(lambda_d0x)*(1-rand(1,sum(I2d)));
    lambda_d0y(I2d) = lambda_d0y(I2d) + shift*max(lambda_d0x)*(1-rand(1,sum(I2d)));
    lambda_Cx(I2C) = lambda_Cx(I2C) + shift*max(lambda_d0x)*(1-rand(1,sum(I2C)));
    lambda_Cy(I2C) = lambda_Cy(I2C) + shift*max(lambda_d0x)*(1-rand(1,sum(I2C)));
    
    n1ldx = bsxfun(@times,n1,lambda_d0x);
    n2ldy = bsxfun(@times,n2,lambda_d0y);

    Lamdx = repmat(lambda_d0x(:),1,m2); Lamdx = Lamdx(:)';
    Lamdy = repmat(lambda_d0y(:),1,m2); Lamdy = Lamdy(:)'; % d0 changes in each sub block
    LamCx = repmat(lambda_Cx(:)',m1,1); LamCx = LamCx(:)';
    LamCy = repmat(lambda_Cy(:)',m1,1); LamCy = LamCy(:)';

    n1Ldx = bsxfun(@times,n1,Lamdx);
    n2Ldy = bsxfun(@times,n2,Lamdy);
    n1LCx = bsxfun(@times,n1,LamCx);
    n2LCy = bsxfun(@times,n2,LamCy);

    Q1 = abs(n1Ldx-n1LCx+n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % - + +
    Q2 = abs(n1Ldx-n1LCx-n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % - - -
    Q3 = abs(n1Ldx-n1LCx+n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % - + -
    Q4 = abs(n1Ldx+n1LCx-n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % + - - 
    Q5 = abs(n1Ldx+n1LCx+n2Ldy-n2LCy)< eps^(2/3)*max(lambda_Cx); % + + -
    Q6 = abs(n1Ldx-n1LCx-n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % - - +
    Q7 = abs(n1Ldx+n1LCx-n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % + - +
    Q8 = abs(n1Ldx+n1LCx+n2Ldy+n2LCy)< eps^(2/3)*max(lambda_Cx); % + + +

    Q1d = abs(n1ldx-n2ldy) < eps^(2/3)*max(lambda_d0x);
    Q2d = abs(n1ldx+n2ldy) < eps^(2/3)*max(lambda_d0x);
    
    I = logical(sum(Q1d+Q2d));
    I2 = logical(sum(Q1+Q2+Q3+Q4+Q5+Q6+Q7+Q8));
    shift = shift*sqrt(10);

    if shift > 0.1*maxVal
        warning('Needed to shfit basis functions by 10% of max')
    end
end

SLambdaC = sig_fC^2*(2*pi)*lCx*lCy*exp(-0.5*(lambda_Cx.^2*lCx^2+lambda_Cy.^2*lCy^2));
SLambdad0 = sig_fd0^2*(2*pi)*ldx*ldy*exp(-0.5*(lambda_d0x.^2*ldx^2+lambda_d0y.^2*ldy^2));

end

