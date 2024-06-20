function [maxA,varargout] = maxNDgrid(A)
% [maxA,ix,iy,iz....] = maxNDgrid(A)
%maxNDGRID This script finds the (ix,iy,iz) location and value (maxA) of
%the maximum of a MxNxPx... matrix, e.g. an error grid where:
% ix is the first dimension index of the maximum - the (ix,...) cooridinate
% iy is the second dimension index of the maximum - the (~,iy,...) coordinate
% iz is the third dimension index of the maximum - the (~,~,iz,...) coordinate
% etc.
% such that maxA = A(ix,iy,iz,...)
Nd = ndims(A);
maxA = max(A,[],'all');
if Nd == 2
    [~,ix,iy] = maxgrid(A);
    % outputs
    varargout = {ix,iy};
    % confirm 
    A(ix,iy)==maxA;
elseif Nd == 3
    % zth dimension
    x = squeeze(max(squeeze(max(A))));
    [~,iz] = max(x);
    % yth dimension
    x = squeeze(max(A(:,:,iz)));
    [~,iy] = max(x);
    % xth dimension
    [~,ix] = max(A(:,iy,iz));
   % outputs
    varargout = {ix,iy,iz};
    % confirm 
    A(ix,iy,iz)==maxA;
end