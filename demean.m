function [ dV ] = demean( V )
%[ dV ] = demean( V )
% quick function to remove mean from columns of V

dV = detrend(V,'constant');


end

