function x = axcorn(ax,icorner,xory)
%  x = axposc(ax,val) 
%   Function to grab corner locations [x,y] within the figure window
%   for a particular set of axes (or corner of the figure on the screen). 
%   So 'x' is by default an Nx2 matrix of x,y corner locations for N
%   corners.
% 
%  This is to allow one-line operations on the axes location.
% 
%   The order of the corners is anticlockwise:
%     1) Bottom Left, 2) Bottom Right, 3) Top Right, 4) Top Left
%   If 'icorner' is blank, give all corners, in that order. 
% 
%  Can request just 'x' coordinates of corner, or just 'y' coordinates of
%  corner, or 'xy' (default) both coordinates of corner.

if nargin < 2 || isempty(icorner)
    icorner = 1:4;
end

if nargin < 3 || isempty(xory)
    xory = 'xy';
end

% seek x or y or both
clm = find(any(xory == ['x';'y'],2)); 

% find corner locations
pos = get(ax,'pos');

cornx = [pos(1),pos(1)+pos(3),pos(1)+pos(3),pos(1)]';
corny = [pos(2),pos(2),pos(2)+pos(4),pos(2)+pos(4)]';
corn = [cornx,corny];

x = corn(icorner,clm); %#ok<FNDSB> 

end