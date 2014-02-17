function [ answ ] = Dang( a1,a2,acrit,unit,axial )
%[ answ ] = Dang( a1,a2,acrit,unit,axial )
% Function to assess whether two angles differ by more or less than a
% critical amount
% INPUTS
% a1    : angle 1
% a2    : angle 2
% acrit : if a1 and a2 differ by more than this, the fn will output 0
% unit  : 'deg' or 'rad' - 
% axial : 'axis' or 'direct' - option to care about polarisation (i.e. are 0 and 180 the same?)

% setup
if nargin < 3
    error('Not enough inputs - need two angles and a critical angle')
end
if nargin < 4
    unit='deg';
end
if nargin < 5
    axial='direct';
end
% et all in right form
if strcmp(unit,'rad')==1
    a1=r2d(a1);
    a2=r2d(a2);
    acrit=r2d(acrit);
end
if strcmp(axial,'axis')==1 % move into 0-180 range and then stretch
    a1=mod(a1,180); a1=2*a1;
    a2=mod(a2,180); a2=2*a2;
    acrit=2*acrit;
end
% processing
da = abs(a1-a2);
mrlR = 0.5*sqrt((1+cosd(da)).^2 + sind(da).^2);
mrlRcrit = 0.5*sqrt((1+cosd(acrit))^2 + sind(acrit)^2);
% output
if mrlR < mrlRcrit 
answ=0;
else
answ=1;
end




end

