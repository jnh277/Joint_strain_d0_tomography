function [ y ] = make_measurements( entry, exit, nhat, Fxx, Fxy, Fyy)
%[ y ] = make_measurements( entry, exit, nhat, Fxx, Fxy, Fyy)
% inputs:
%   entry: the ith colum contains [x, y]' for the entry point of ray i
%   exit: the ith column contains [x,y]' for the exit point of the ith ray
%   nhat: the ith column contains the unit vector direction for the ith ray
%   Fxx: an anonymous function for the exx strain field Exx = Fxx(x,y)
%   Fxy: an anonymous function for the exx strain field Exy = Fxy(x,y)
%   Fyy: an anonymous function for the exx strain field Eyy = Fxy(x,y)

i_res = 500;
% create a whole bunch of points that go from entry to exit
Rx = linspace2D(entry(1,:)',exit(1,:)',i_res );
Ry = linspace2D(entry(2,:)',exit(2,:)',i_res );
% evaluate the strain fields at these points
Vxx = Fxx(Rx,Ry);
Vxy = Fxy(Rx,Ry);
Vyy = Fyy(Rx,Ry);
N1 = repmat(nhat(1,:)'.^2,1,i_res);
N12 = repmat(2*nhat(1,:)'.*nhat(2,:)',1,i_res);
N2 = repmat(nhat(2,:)'.^2,1,i_res);
% do math
V =  N1.*Vxx+N12.*Vxy+N2.*Vyy;
% S = sqrt((Rx-Rx(:,1)).^2 + (Ry-Ry(:,1)).^2);
% L = S(:,end) - S(:,1);
% ds = L/i_res;
% integrate
y = trapz(1:i_res,V,2)/(i_res-1);

end

