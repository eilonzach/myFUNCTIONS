function [ dV,Vav ] = demean( V )
%[ dV , Vav] = demean( V )
% quick function to remove mean from columns of V

Vav = nanmean(V);

inan = isnan(V);
dV(~inan) = detrend(V(~inan),'constant');
dV(inan) = nan;

end

