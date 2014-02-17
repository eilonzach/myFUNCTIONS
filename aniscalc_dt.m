function [ dT ] = aniscalc_dt( anis, Vav, L )
%[ dT ] = aniscalc_dt( anis, Vav, L )
%   
% Function to calculate the expected time delay between fast and slow
% quasi-wavelets that have passed through an anisotropic layer of thickness
% "L", with average velocity "Vav" and anisotropy "anis"
%
% INPUTS:
% anis      - the anisotropy of the layer, given as a percentage (e.g. for
%               a 7% anisotropic medium, anis=7
% Vav       - average shear wave velocity of the layer (km/s)
% L         - thickenss of the layer (km)
%
% EQUATIONS: (for reference)
% 
% (1)  Vav  = (v1 + v2)/2        [ v1 is fast, v2 is slow ]
% (2)  anis = (v1 - v2)/Vav
% (3)  dT   = L/v2 - L/v1
%
% from (1) and (2) we also have
% (4)  v1   = Vav*(1+anis)
% (5)  v2   = Vav*(1-anis)

% By simply combining (3),(4),(5) above, we obtain:

A = anis./100;

dT = (L./Vav).*A./(1 - 0.25.*A.^2);

% Thus clearly L = dT*(Vav./A)*(4 - A.^2)/4;


end

