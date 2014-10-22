function [rgbs] = colour_get(v,vmax,vmin,colourmap)
%  [rgbs] = colour_get(v,vmax,vmin,colourscheme)
% 
% This function takes a vector of numerical values and returns the
% corresponding [r g b] values so they can be plotted by colour
%
% INPUT	  v 	Nx1 vector of numerical values
% 		  vmax	maximum on scale (optional, will use 1 as default)
% 		  vmin	minimum on scale (optional, will use 0 as default)
%		  colourscheme 	defines name of colormap from built-in options (optional) 
% N.B. values of v outside the defined scale will be returned as the colours of the extremes
%
% OUTPUT  rgbs	Nx3 matrix specifying colours corresponding to numerical values in v
% N.B. uses colormap to get vectors and then linear interpolates to values

if nargin < 2
vmax = max(v);
vmin = min(v);
end

if nargin < 4
cmap = colormap;
else
cmap = colourmap;
end

v(v>vmax) = vmax;
v(v<vmin) = vmin;

M = length(cmap);
N = length(v);
Dv = vmax-vmin; 
if Dv <= 0, error('vmax must be larger than vmin'); end

rgbs = zeros(N,3);
cvalues = vmin:Dv/(M-1):vmax;
for ic = 1:3 % loop over RGB
rgbs(:,ic) = interp1(cvalues,cmap(:,ic),v);
end 

end