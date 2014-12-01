function [ indxs ] = balanced_resampling( Ndata,niter )
%[ indxs ] = balanced_resampling( Ndata,niter )
%   function to get indices for balanced resampling of data for
%   bootstrapping
%
% INPUTS:
%   Ndata   = number of data in original dataset
%   niter   = number of iterations in bootstrap
%   
% OUTPUTS:
%   indxs   = Ndata x niter matrix, each column of which is a set of
%              indices from the original dataset. Each index appears a
%              total of precisely niter times in the whole matrix

xx = [1:Ndata]';

XX = repmat(xx,niter,1);

ind = randperm(Ndata*niter);

YY = XX(ind);

indxs = reshape(YY,Ndata,niter);



end

