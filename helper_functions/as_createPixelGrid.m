function [ pixelPoints ] = as_createPixelGrid(hRes,width)
% [pixelPoints] = as_createPixelGrid(hRes,width)
%
% USAGE:
% General
%
% INPUTS:
% hRes - number of detector pixels
% width - detector width
%
% OUTPUTS:
% pixelPoints - (x,y) for each detector pixel
%
% PURPOSE:
% Code Description:
% Generates detector pixel co-ordinates for a given resolution and width.
% Alexander Gregg, The University of Newcastle, 2016
%----------------------------------------------------------------------

% Calculate the distance between pixels.
hDiv = width/(hRes);
% Calculate the pixel center co-ordinates.
xgv = (-width/2+hDiv/2):hDiv:(width/2-hDiv/2);
% Return the pixel center co-ordinates.
pixelPoints =  [zeros(length(xgv),1),xgv']; 
end

