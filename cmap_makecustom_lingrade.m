function [ cc ] = cmap_makecustom_lingrade( c1, c2,varargin )
%[ cc ] = cmap_makecustom_lingrade( c1, c2, .... )
%   Makes custom colourmap scaling linearly between as many colours as are
%   input (must be at least two)
% 
% INPUT:
%   c1      First endmember colour as 3-component vector [r g b]
%   c2      Second endmember colour as 3-component vector [r g b]
%   ...
% 
% OUTPUT:
%   cc      colourmap: a (200 x 3) matrix, can be called with colormap(cc)
% 
% Written by Zach Eilon, 2020

cc_ends = [c1(:)';c2(:)'];
if nargin >= 3
    for ic = 1:length(varargin)
        cc_ends = [cc_ends;varargin{ic}(:)'];
    end
end

N = 201; % nominally 200 long (may actually be a bit longer if no simple division; let matlab do the interpolation later)
Nc = size(cc_ends,1);

M = ceil((N-1)/(Nc-1)); 

% prep colour matrix
cc = zeros(M*(Nc-1) + 1,3);

rampu = linspace(0,1,M+1)'; rampu = rampu(1:end-1);
rampd = linspace(1,0,M+1)'; rampd = rampd(1:end-1);

for ic = 1:Nc-1
    j1 = 1 + (ic-1)*M;
    j2 = ic*M;
    
    cc(j1:j2,:) = rampd*cc_ends(ic,:) + rampu*cc_ends(ic+1,:);
end
cc(end,:) = cc_ends(end,:);


end

