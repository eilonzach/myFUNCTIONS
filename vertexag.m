function [ Vfac ] = vertexag( ax )
%[ Vfac ] = vertexag( ax )
%
% function to compute the vertical exaggeration by empirically comparing
% the vertical pixels (Y) vs. horizontal pixels (X) to 
% the vertical scale  (y) vs. horizontal scale  (x).

if nargin < 1 || isempty(ax)
    ax = gca;
end

x = diff(get(ax,'xlim'));
y = diff(get(ax,'ylim'));
bx = get(ax,'PlotBoxAspectRatio');
X = bx(1);
Y = bx(2);

Vfac = (x/X)./(y/Y);



end

