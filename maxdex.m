function [ maxX_ind ] = maxdex( X,a )
% [ maxX_ind ] = maxdex( X,a )
%   simple function to return the index of the maximum point in vector X
%   
%   if a second argument is given, the function outputs the index of the
%   point in X furthest from the water level, "a"
%   basically just outputs the second output of the "max" function,
%   without giving you the magnitude of the maximum value.
% 
%   N.B. can be used as a zero finder if a==0
%
%   Intended for use when calling the value in one vector corresponding to
%   the maximum value in X - i.e more efficient than the clunkier:
%       Y(find(X==max(X)) or, more often, Y(find((X-a)==max(X-a)))
%   Instead, can now use
%       Y(maxdex(X)) or Y(maxdex(X,a))
%   See also mindex.m

if nargin<2
    [~,maxX_ind] = max(X);
else
    [~,maxX_ind] = max(abs(X-a));
end

end