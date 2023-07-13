function [ xx ] = fix_level( x, y )
%[xx] = fix_level(x,level)
% fixes the number/vector x to the nearest y (the level)
xx = fix(x./y).*y;
end