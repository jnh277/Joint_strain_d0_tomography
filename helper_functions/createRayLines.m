function [ lines ] = createRayLines(angles, pixelPoints)
% [ lines ] = as_createRayLines(angles, pixelPoints,numAngles,numPixels)
%
% USAGE:
% General
%
% INPUTS:
% angles - vector of projection angles
% pixelPoints - (x,y) for each detector pixel
%
% OUTPUTS:
% lines - ray origin and direction co-ordinates (pixel, pixel + nhat)
%
% PURPOSE:
% Code Description:
% Generates geometry for each unique neutron ray.
% Alexander Gregg, The University of Newcastle, 2016
% Johannes Hendriks, The University of Newcastle, 2017
%----------------------------------------------------------------------

nHat = [1;0;0];
[numPixels,~] = size(pixelPoints);
[numAngles] = numel(angles);

% Preallocate the lines hypermatrix for speed. 
lines = NaN(2,2,numPixels*numAngles);
% Loop over each angle
for i = 1:numAngles
    % Create a rotation matrix using the ray propagation angle.
    R = [cos(angles(i)), -sin(angles(i)), 0;
        sin(angles(i)), cos(angles(i)), 0
        0, 0, 1];  
%     rotationMat([0,0, angles(i)]);
    % Extract the components pertinent to 2D.
    R2d = R(1:2,1:2);
    % Loop over each pixel within each angle.
   for j = 1:numPixels
       % Assign a unique index to this angle-pixel combination.
       ind = (i-1)*numPixels + j;
       % Rotate the ray propagation vector to the new angle.
       RnHat = R*nHat;
       % Create the line object.
       lines(:,:,ind) = [R2d*(pixelPoints(j,:)'), R2d*(pixelPoints(j,:)')+RnHat(1:2)];    
   end
end
end


