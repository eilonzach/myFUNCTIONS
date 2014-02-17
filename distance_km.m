function [ KMLEN, AZ ] = distance_km(LAT1,LON1,LAT2,LON2)
%[ KMLEN, AZ ] = distance_km(LAT1,LON1,LAT2,LON2)
%   This is just the distance function but with the Earth's ellipsoid input
%   by defult, such that the distance output is in km.

[KMLEN,AZ] = distance(LAT1,LON1,LAT2,LON2,[6378.137, 0.08181919]);



end

