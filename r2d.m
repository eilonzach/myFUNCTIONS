function [ deg ] = r2d( rad )
% Converts input in radians to output in degrees
if nargin==0
    rad=1;
end
deg=rad*180/pi;
end