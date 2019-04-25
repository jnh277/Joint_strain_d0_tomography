function [ d ] = make_measurements_d0varying( entry, exit,Fxx, Fxy, Fyy,d0_field)
%[ y ] = make_measurements( entry, exit, nhat, Fxx, Fxy, Fyy)
% inputs:
%   entry: the ith colum contains [x, y]' for the entry point of ray i
%   exit: the ith column contains [x,y]' for the exit point of the ith ray
%   nhat: the ith column contains the unit vector direction for the ith ray
%   Fxx: an anonymous function for the exx strain field Exx = Fxx(x,y)
%   Fxy: an anonymous function for the exx strain field Exy = Fxy(x,y)
%   Fyy: an anonymous function for the exx strain field Eyy = Fxy(x,y)
%   d0_field: an anonymous function for the d0 field d0 = d0_field(x,y)

n = exit - entry;
L = hypot(n(1,:),n(2,:));
nhat = n./L;

N1 = nhat(1,:).^2;
N12 = 2*nhat(1,:).*nhat(2,:);
N2 = nhat(2,:).^2;
% x = entry(1,:)+n(1,:)*tau
integrand = @(tau) (N1.*Fxx(entry(1,:)+n(1,:)*tau,entry(2,:)+n(2,:)*tau)+...
    N12.*Fxy(entry(1,:)+n(1,:)*tau,entry(2,:)+n(2,:)*tau)...
    +N2.*Fyy(entry(1,:)+n(1,:)*tau,entry(2,:)+n(2,:)*tau)).*d0_field(entry(1,:)+n(1,:)*tau,entry(2,:)+n(2,:)*tau)...
    + d0_field(entry(1,:)+n(1,:)*tau,entry(2,:)+n(2,:)*tau);

% integrand2 = @(tau) (N1(1).*Fxx(entry(1,1)+nhat(1,1)*tau,entry(2,1)+nhat(2,1)*tau)+...
%     N12(1).*Fxy(entry(1,1)+nhat(1,1)*tau,entry(2,1)+nhat(2,1)*tau)...
%     +N2(1).*Fyy(entry(1,1)+nhat(1,1)*tau,entry(2,1)+nhat(2,1)*tau)).*d0_field(entry(1,1)+nhat(1,1)*tau,entry(2,1)+nhat(2,1)*tau)...
%     + d0_field(entry(1,1)+nhat(1,1)*tau,entry(2,1)+nhat(2,1)*tau);

% d2 = integral(integrand2,0,L(1))/L(1)

d = integral(integrand,0,1,'ArrayValued',true)';        % .*L for the 0 to 1 scaling and ./L for the average


end

