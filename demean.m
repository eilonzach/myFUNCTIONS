function [ V_,Vav ] = demean( V,opt )
%[ dV , Vav] = demean( V,opt )
% quick function to remove mean from columns of V
% "opt" is 'constant' (default) or the N of the polynomial of the trend to 
% fit and remove

if nargin < 2 || isempty(opt)
    opt = 'constant'; % default is just to remove a constant value
end

Vav = nanmean(V);

V_ = nan(size(V));

% inan = isnan(V);
try
    V_ = detrend(V,opt,'omitnan');
catch
    warning('detrend() with nans breaking for some reason');
    V_col_nnan = ~isnan(nanmean(V));
    V_(:,V_col_nnan) = detrend(V(:,V_col_nnan),opt,'omitnan');
end
% V_(inan) = nan;

end

