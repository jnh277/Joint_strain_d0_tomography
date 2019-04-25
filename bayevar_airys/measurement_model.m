function [ypred] = measurement_model(ipbf,id0bf,w_d0,w_C,nhat,L,nu)
%ypred] = predict_y(ipbf,w_d0,w_C,nhat,L,nu)
%   where we will vectorise [w_{d0,k}*w_{C,j}] assuming that k is the row
%   index

W1W2 = w_d0(:).*w_C(:)';
W1W2 = W1W2(:);

n1 = nhat(1,:)';
n2 = nhat(2,:)';
n1n1 = n1.^2;
n1n2 = 2*n1.*n2;
n2n2 = n2.^2;

idxdxV = ipbf(:,:,1)*W1W2;
idydyV = ipbf(:,:,2)*W1W2;
idxdyV = ipbf(:,:,3)*W1W2;

ypred = (n1n1.*(idydyV - nu*idxdxV) -(1+nu)*n1n2.*idxdyV + n2n2.*(idxdxV-nu*idydyV) + id0bf*w_d0)./L(:);

end

