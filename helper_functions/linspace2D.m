function [ r ] = linspace2D(x0,x1,N )
%[ r ] = linspace2D(x0,x1,N )
%   efficiently does the linspace between two column vectors

assert(iscolumn(x0) && iscolumn(x1))

dx = (x1-x0)/(N-1);
r = repmat(dx,1,N);
r(:,1) = x0;
r = cumsum(r,2);


end

