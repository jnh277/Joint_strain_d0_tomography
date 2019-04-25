function [ lines ] = gen_rays_3d(angles, num_pixels, width)
% [ lines ] = as_createRayLines(angles, pixelPoints,numAngles,numPixels)
%
% USAGE:
% General
%
% INPUTS:
% angles - 3XN array containing the euler angles for each rotation [phi,
% theta, psi] or [roll,pitch,yaw], or [about x,about y, about z]
% num_pixels - number of pixels in one direction (assumed grid is square)
% width - width of the grid (assumed square)
%
% OUTPUTS:
% lines - ray origin and direction co-ordinates ([pixel, nhat]')
%
%----------------------------------------------------------------------

hDiv = width/(num_pixels);

% Calculate the pixel center co-ordinates.
xgv = (-width/2+hDiv/2):hDiv:(width/2-hDiv/2);

[y_grid,z_grid] = meshgrid(xgv);
grid_points = [zeros(1,num_pixels^2);y_grid(:)';z_grid(:)'];

nHat = [1;0;0];

[~,numAngles] = size(angles);

% Preallocate the lines hypermatrix for speed. 
lines = NaN(6,num_pixels^2*numAngles);
% Loop over each angle
for i = 1:numAngles
    % Create a rotation matrix using the ray propagation angle.
    R = eulerRotation(angles(:,i));
    lines(1:3,(i-1)*num_pixels^2+1:i*num_pixels^2) = R*grid_points;
    lines(4:6,(i-1)*num_pixels^2+1:i*num_pixels^2) = repmat(R*nHat,1,num_pixels^2);
end
end


