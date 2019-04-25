function [ipbf_order1,ipbf_order2] = intprodbasisfunc(lambda_d0x,lambda_d0y,L_d0x,L_d0y,...
    lambda_Cx,lambda_Cy,L_Cx,L_Cy,entry,exit,nhat)
%UNTITLED8 Summary of this function goes here

% order refers to whether when vectorising [w_{d0,k}*w_{C,j}] we take k or
% j as the row index,
% if order == 1 then we use k as the row index, which means d0 basis
% changes withing sub blocks
% if order == 2 then we use j as the row index, which means C basis changes
% within sub blocks

%   Detailed explanation goes here
    m1 = length(lambda_d0x);
    m2 = length(lambda_Cx);
    n = size(entry,2);
    
    Lamdx = repmat(lambda_d0x(:),1,m2); Lamdx = Lamdx(:)';
    Lamdy = repmat(lambda_d0y(:),1,m2); Lamdy = Lamdy(:)'; 
    LamCx = repmat(lambda_Cx(:)',m1,1); LamCx = LamCx(:)';
    LamCy = repmat(lambda_Cy(:)',m1,1); LamCy = LamCy(:)';

    % at entry
    alpha_Cx = LamCx.*(entry(1,:)'+L_Cx);
    alpha_Cy = LamCy.*(entry(2,:)'+L_Cy);
    alpha_d0x = Lamdx.*(entry(1,:)'+L_d0x);
    alpha_d0y = Lamdy.*(entry(2,:)'+L_d0y);
    
    Gamma10 = sin(alpha_d0x-alpha_Cx+alpha_d0y+alpha_Cy); % - + +
    Gamma20 = sin(alpha_d0x-alpha_Cx-alpha_d0y-alpha_Cy); % - - -
    Gamma30 = sin(alpha_d0x-alpha_Cx+alpha_d0y-alpha_Cy); % - + -
    Gamma40 = sin(alpha_d0x+alpha_Cx-alpha_d0y-alpha_Cy); % + - - 
    Gamma50 = sin(alpha_d0x+alpha_Cx+alpha_d0y-alpha_Cy); % + + -
    Gamma60 = sin(alpha_d0x-alpha_Cx-alpha_d0y+alpha_Cy); % - - +
    Gamma70 = sin(alpha_d0x+alpha_Cx-alpha_d0y+alpha_Cy); % + - +
    Gamma80 = sin(alpha_d0x+alpha_Cx+alpha_d0y+alpha_Cy); % + + +
    
    % at exit
    alpha_Cx = LamCx.*(exit(1,:)'+L_Cx);
    alpha_Cy = LamCy.*(exit(2,:)'+L_Cy);
    alpha_d0x = Lamdx.*(exit(1,:)'+L_d0x);
    alpha_d0y = Lamdy.*(exit(2,:)'+L_d0y);
    
    Gamma11 = sin(alpha_d0x-alpha_Cx+alpha_d0y+alpha_Cy); % - + +
    Gamma21 = sin(alpha_d0x-alpha_Cx-alpha_d0y-alpha_Cy); % - - -
    Gamma31 = sin(alpha_d0x-alpha_Cx+alpha_d0y-alpha_Cy); % - + -
    Gamma41 = sin(alpha_d0x+alpha_Cx-alpha_d0y-alpha_Cy); % + - - 
    Gamma51 = sin(alpha_d0x+alpha_Cx+alpha_d0y-alpha_Cy); % + + -
    Gamma61 = sin(alpha_d0x-alpha_Cx-alpha_d0y+alpha_Cy); % - - +
    Gamma71 = sin(alpha_d0x+alpha_Cx-alpha_d0y+alpha_Cy); % + - +
    Gamma81 = sin(alpha_d0x+alpha_Cx+alpha_d0y+alpha_Cy); % + + +
    
    n1 = nhat(1,:)';
    n2 = nhat(2,:)';
    
    n1Ldx = bsxfun(@times,n1,Lamdx);
    n2Ldy = bsxfun(@times,n2,Lamdy);
    n1LCx = bsxfun(@times,n1,LamCx);
    n2LCy = bsxfun(@times,n2,LamCy);
    
    Gamma1 = (Gamma11 - Gamma10)./(8*(n1Ldx-n1LCx+n2Ldy+n2LCy)); % - + +
    Gamma2 = (Gamma21 - Gamma20)./(8*(n1Ldx-n1LCx-n2Ldy-n2LCy)); % - - -
    Gamma3 = (Gamma31 - Gamma30)./(8*(n1Ldx-n1LCx+n2Ldy-n2LCy)); % - + -
    Gamma4 = (Gamma41 - Gamma40)./(8*(n1Ldx+n1LCx-n2Ldy-n2LCy)); % + - - 
    Gamma5 = (Gamma51 - Gamma50)./(8*(n1Ldx+n1LCx+n2Ldy-n2LCy)); % + + -
    Gamma6 = (Gamma61 - Gamma60)./(8*(n1Ldx-n1LCx-n2Ldy+n2LCy)); % - - +
    Gamma7 = (Gamma71 - Gamma70)./(8*(n1Ldx+n1LCx-n2Ldy+n2LCy)); % + - +
    Gamma8 = (Gamma81 - Gamma80)./(8*(n1Ldx+n1LCx+n2Ldy+n2LCy)); % + + +
    
    Omega1 = (-Gamma1-Gamma2+Gamma3+Gamma4-Gamma5+Gamma6-Gamma7+Gamma8);
    Omega2 = (-Gamma1+Gamma2-Gamma3+Gamma4-Gamma5+Gamma6+Gamma7-Gamma8);
    
    ipbf_order1 = nan(n,m1*m2,3);
    % integral of (d/dxdx C) * d0
    ipbf_order1(:,:,1) = (-LamCx.^2/sqrt(L_Cx*L_Cy*L_d0x*L_d0y)).*Omega1;
    % integral of (d/dydy C) * d0
    ipbf_order1(:,:,2) = (-LamCy.^2/sqrt(L_Cx*L_Cy*L_d0x*L_d0y)).*Omega1;
    % integral of (d/dxdy C) * d0
    ipbf_order1(:,:,3) = ((LamCy.*LamCx)/sqrt(L_Cx*L_Cy*L_d0x*L_d0y)).*Omega2;
    
    if nargout == 2
        ipbf_order2 = nan(n,m1*m2,3);
        T1 = reshape(ipbf_order1(:,:,1),n,m1,m2);
        T2 = permute(T1,[1,3,2]);
        ipbf_order2(:,:,1) = reshape(T2,n,m1*m2);
        
        T1 = reshape(ipbf_order1(:,:,2),n,m1,m2);
        T2 = permute(T1,[1,3,2]);
        ipbf_order2(:,:,2) = reshape(T2,n,m1*m2);
        
        T1 = reshape(ipbf_order1(:,:,3),n,m1,m2);
        T2 = permute(T1,[1,3,2]);
        ipbf_order2(:,:,3) = reshape(T2,n,m1*m2);
    end
end

