function [ xx ] = ceil_level( x, y )
%[xx] = ceil_level(x,level)
% ceils the number/vector x to the nearest y (the level)
xx = ceil(x./y).*y;
end

