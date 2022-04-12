function tapr = cos_taper_hipass(k_in,hipass_km,lostop_km)
% tapr = cos_taper_hipass(k_in,hipass_km,lostop_km)
%
% function to create a cosine taper in wavenumber space. The inputs are a
% matrix of wavenumbers (k_in, in units of 2pi/km), a high-pass wavelength
% (hipass_km, in units of km), and a low-stop wavelength (lostop_km, in
% units of km). The filter acts to completely retain wavelengths shorter
% (=wavenumbers larger) than hipass_km, and completely remove wavelengths
% longer (=wavenumbers smaller) than lostop_km. Between these two
% wavelengths, there is a cosine taper in wavenumber space that grades from
% 1 to 0. The output, tapr, is a matrix the same size as k_in, with
% numbers from 1 (retained) to 0 (blocked) corresponding to the wavenumber
% values in k_in.

if hipass_km > lostop_km
    error('the high pass must be shorter wavelength than the lostop')
end

% relevant wavenumber increment to think about
dk = mean(mean(abs(diff(k_in))));
% convert filter edges to wavenumbners from wavelengths
hipass_k = 2*pi./hipass_km;
lostop_k = 2*pi./lostop_km;
% make index vector of wavenumbers
k_index = unique([min(min(k_in)):dk:hipass_k,hipass_k:dk/10:lostop_k,lostop_k:dk:max(max(k_in)),max(max(k_in))]');
% make corresponding index vector of taper (1 - 0)
tapr_index = ones(size(k_index));
tapr_index(k_index<lostop_k) = 0;
w = tukeywin(101,1); w = w(1:51);
tapr_index( k_index >= lostop_k & k_index <= hipass_k ) ...
    = interp1(linspace(lostop_k,hipass_k,51),w,k_index( k_index >= lostop_k & k_index <= hipass_k ) );

% make actual taper for each k_in by interpolating with the index vectors
tapr = interp1(k_index,tapr_index,k_in,'linear',0);

end
