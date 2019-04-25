function [Jwd0,JwC] = measurement_model_J(ipbf_order1,ipbf_order2,id0bf,w_d0,w_C,nhat,L,nu)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here

n1 = nhat(1,:)';
n2 = nhat(2,:)';
n1n1L = n1.^2./L(:);
n1n2L = 2*n1.*n2./L(:);
n2n2L = n2.^2./L(:);

m1 = length(w_d0);
m2 = length(w_C);
n = size(nhat,2);       % number of measurements

% for derivatives with respect to w_C using order 1
Wd0 = w_d0(:).*ones(1,m2);
Wd0 = Wd0(:)';

Mxx = ipbf_order1(:,:,1).*Wd0;       % using singleton expansion and then a block sum
Myy = ipbf_order1(:,:,2).*Wd0;
Mxy = ipbf_order1(:,:,3).*Wd0;

mxx = blockproc(Mxx,[n,m1],@(x)sum(x.data,2));
myy = blockproc(Myy,[n,m1],@(x)sum(x.data,2));
mxy = blockproc(Mxy,[n,m1],@(x)sum(x.data,2));

JwC = n1n1L.*(myy -nu*mxx)-(1+nu)*n1n2L.*mxy+n2n2L.*(mxx-nu*myy);        % the 1/L was done earlier


% for derivatives with respect to w_d0 we need to use order 2 (helps with
% blockpro)

WC = w_C(:).*ones(1,m1);
WC = WC(:)';

Mxx = ipbf_order2(:,:,1).*WC;       % using singleton expansion and then a block sum
Myy = ipbf_order2(:,:,2).*WC;
Mxy = ipbf_order2(:,:,3).*WC;

mxx = blockproc(Mxx,[n,m2],@(x)sum(x.data,2));
myy = blockproc(Myy,[n,m2],@(x)sum(x.data,2));
mxy = blockproc(Mxy,[n,m2],@(x)sum(x.data,2));

Jwd0 = n1n1L.*(myy -nu*mxx)-(1+nu)*n1n2L.*mxy+n2n2L.*(mxx-nu*myy) + id0bf./L(:); 
end

