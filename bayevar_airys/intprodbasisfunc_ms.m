function [ipbf_order1,ipbf_order2] = intprodbasisfunc_ms(lambda_d0x,lambda_d0y,L_d0x,L_d0y,...
    lambda_Cx,lambda_Cy,L_Cx,L_Cy,entry,exit,nhat,nsegs)
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

n1 = nhat(1,:)';
n2 = nhat(2,:)';

n1Ldx = bsxfun(@times,n1,Lamdx);
n2Ldy = bsxfun(@times,n2,Lamdy);
n1LCx = bsxfun(@times,n1,LamCx);
n2LCy = bsxfun(@times,n2,LamCy);

Gamma1_seg = zeros(n,m1*m2);
Gamma2_seg = zeros(n,m1*m2);
Gamma3_seg = zeros(n,m1*m2);
Gamma4_seg = zeros(n,m1*m2);
Gamma5_seg = zeros(n,m1*m2);
Gamma6_seg = zeros(n,m1*m2);
Gamma7_seg = zeros(n,m1*m2);
Gamma8_seg = zeros(n,m1*m2);

for ss = 1:max(nsegs)
    segInd = find(nsegs >= ss);        % index of measurements that have at least this many segments
    
    x0 = entry(ss*2-1,segInd)'; y0 = entry(ss*2,segInd)';  
    % at entry
    alpha_Cx = LamCx.*(x0+L_Cx);
    alpha_Cy = LamCy.*(y0+L_Cy);
    alpha_d0x = Lamdx.*(x0+L_d0x);
    alpha_d0y = Lamdy.*(y0+L_d0y);
    
    Gamma10 = sin(alpha_d0x-alpha_Cx+alpha_d0y+alpha_Cy); % - + +
    Gamma20 = sin(alpha_d0x-alpha_Cx-alpha_d0y-alpha_Cy); % - - -
    Gamma30 = sin(alpha_d0x-alpha_Cx+alpha_d0y-alpha_Cy); % - + -
    Gamma40 = sin(alpha_d0x+alpha_Cx-alpha_d0y-alpha_Cy); % + - - 
    Gamma50 = sin(alpha_d0x+alpha_Cx+alpha_d0y-alpha_Cy); % + + -
    Gamma60 = sin(alpha_d0x-alpha_Cx-alpha_d0y+alpha_Cy); % - - +
    Gamma70 = sin(alpha_d0x+alpha_Cx-alpha_d0y+alpha_Cy); % + - +
    Gamma80 = sin(alpha_d0x+alpha_Cx+alpha_d0y+alpha_Cy); % + + +
    
    % at exit
    xf = exit(ss*2-1,segInd)'; yf = exit(ss*2,segInd)';  
    alpha_Cx = LamCx.*(xf+L_Cx);
    alpha_Cy = LamCy.*(yf+L_Cy);
    alpha_d0x = Lamdx.*(xf+L_d0x);
    alpha_d0y = Lamdy.*(yf+L_d0y);
    
    Gamma11 = sin(alpha_d0x-alpha_Cx+alpha_d0y+alpha_Cy); % - + +
    Gamma21 = sin(alpha_d0x-alpha_Cx-alpha_d0y-alpha_Cy); % - - -
    Gamma31 = sin(alpha_d0x-alpha_Cx+alpha_d0y-alpha_Cy); % - + -
    Gamma41 = sin(alpha_d0x+alpha_Cx-alpha_d0y-alpha_Cy); % + - - 
    Gamma51 = sin(alpha_d0x+alpha_Cx+alpha_d0y-alpha_Cy); % + + -
    Gamma61 = sin(alpha_d0x-alpha_Cx-alpha_d0y+alpha_Cy); % - - +
    Gamma71 = sin(alpha_d0x+alpha_Cx-alpha_d0y+alpha_Cy); % + - +
    Gamma81 = sin(alpha_d0x+alpha_Cx+alpha_d0y+alpha_Cy); % + + +
    
    
    Gamma1_seg(segInd,:) = Gamma1_seg(segInd,:) + (Gamma11 - Gamma10);
    Gamma2_seg(segInd,:) = Gamma2_seg(segInd,:) + (Gamma21 - Gamma20);
    Gamma3_seg(segInd,:) = Gamma3_seg(segInd,:) + (Gamma31 - Gamma30);
    Gamma4_seg(segInd,:) = Gamma4_seg(segInd,:) + (Gamma41 - Gamma40);
    Gamma5_seg(segInd,:) = Gamma5_seg(segInd,:) + (Gamma51 - Gamma50);
    Gamma6_seg(segInd,:) = Gamma6_seg(segInd,:) + (Gamma61 - Gamma60);
    Gamma7_seg(segInd,:) = Gamma7_seg(segInd,:) + (Gamma71 - Gamma70);
    Gamma8_seg(segInd,:) = Gamma8_seg(segInd,:) + (Gamma81 - Gamma80);
    

end
Gamma1 = (Gamma1_seg)./(8*(n1Ldx-n1LCx+n2Ldy+n2LCy)); % - + +
Gamma2 = (Gamma2_seg)./(8*(n1Ldx-n1LCx-n2Ldy-n2LCy)); % - - -
Gamma3 = (Gamma3_seg)./(8*(n1Ldx-n1LCx+n2Ldy-n2LCy)); % - + -
Gamma4 = (Gamma4_seg)./(8*(n1Ldx+n1LCx-n2Ldy-n2LCy)); % + - - 
Gamma5 = (Gamma5_seg)./(8*(n1Ldx+n1LCx+n2Ldy-n2LCy)); % + + -
Gamma6 = (Gamma6_seg)./(8*(n1Ldx-n1LCx-n2Ldy+n2LCy)); % - - +
Gamma7 = (Gamma7_seg)./(8*(n1Ldx+n1LCx-n2Ldy+n2LCy)); % + - +
Gamma8 = (Gamma8_seg)./(8*(n1Ldx+n1LCx+n2Ldy+n2LCy)); % + + +

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

