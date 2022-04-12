function tapr = cos_taper_lopass(k_in,lopass_km,histop_km)
% tapr = cos_taper_lopass(k_in,lopass_km,histop_km)
%
% function to create a cosine taper in wavenumber space. The inputs are a
% matrix of wavenumbers (k_in, in units of 2pi/km), a low-pass wavelength
% (lopass_km, in units of km), and a high-stop wavelength (histop_km, in
% units of km). The filter acts to completely retain wavelengths longer
% (=wavenumbers shorter) than lopass_km, and completely remove wavelengths
% shorter (=wavenumbers larger) than histop_km. Between these two
% wavelengths, there is a cosine taper in wavenumber space that grades from
% 1 to 0. The output, tapr, is a matrix the same size as k_in, with
% numbers from 1 (retained) to 0 (blocked) corresponding to the wavenumber
% values in k_in.

if lopass_km < histop_km
    error('the low pass must be longer wavelength than the histop')
end

% relevant wavenumber increment to think about
dk = max(max(abs(diff(k_in))));
% convert filter edges to wavenumbners from wavelengths
lopass_k = 2*pi./lopass_km;
histop_k = 2*pi./histop_km;
% make index vector of wavenumbers
k_index = unique([min(min(k_in)):dk:lopass_k,lopass_k:dk/10:histop_k,histop_k:dk:max(max(k_in)),max(max(k_in))]');
% make corresponding index vector of taper (1 - 0)
tapr_index = ones(size(k_index));
tapr_index(k_index>histop_k) = 0;
w = tukeywin(101,1); w = w(51:end);
tapr_index( k_index >= lopass_k & k_index <= histop_k) ...
    = interp1(linspace(lopass_k,histop_k,51),w,k_index( k_index >= lopass_k & k_index <= histop_k));

% make actual taper for each k_in by interpolating with the index vectors
tapr = interp1(k_index,tapr_index,k_in,'linear',0);

end
