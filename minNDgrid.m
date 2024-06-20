function [minA,varargout] = minNDgrid(A)
% [minA,ix,iy,iz....] = minNDgrid(A)
%MINNDGRID This script finds the (ix,iy,iz) location and value (minA) of
%the minimum of a MxNxPx... matrix, e.g. an error grid where:
% ix is the first dimension index of the minimum - the (ix,...) cooridinate
% iy is the second dimension index of the minimum - the (~,iy,...) coordinate
% iz is the third dimension index of the minimum - the (~,~,iz,...) coordinate
% etc.
% such that minA = A(ix,iy,iz,...)
Nd = ndims(A);
minA = min(A,[],'all');
if Nd == 2
    [~,ix,iy] = mingrid(A);
    % outputs
    varargout = {ix,iy};
    % confirm 
    A(ix,iy)==minA;
elseif Nd == 3
    % zth dimension
    x = squeeze(min(squeeze(min(A))));
    [~,iz] = min(x);
    % yth dimension
    x = squeeze(min(A(:,:,iz)));
    [~,iy] = min(x);
    % xth dimension
    [~,ix] = min(A(:,iy,iz));
   % outputs
    varargout = {ix,iy,iz};
    % confirm 
    A(ix,iy,iz)==minA;
end