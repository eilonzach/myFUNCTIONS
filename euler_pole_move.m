function [ PP,azs,arclens ] = euler_pole_move( E, V, P, T, dT )
% [ PP ] = euler_pole_move( E, V, P, T, dT )
% function to calculate motion about Euler pole in increments of time
% each increment, the vector of instantaneous motion is calculated and then
% moved along for time dT, before updating to a new location and repeating
% 
% INPUTS
%   E       location of Euler pole, [lat lon]
%   V       rotation rate in degrees/Myr
%   P       starting point, [lat lon]
%   T       total length of time in Myr
%   dT      time increment in Myr
% 
% OUTPUTS
%   PP      matrix of points representing motion [(N+1)*Lat (N+1)*Lon]
%   azs     vector of azimuths between points (Nx1) (deg from north)
%   arclens vector of lengths between points (Nx1)  (degrees)

if nargin<4
    T = 1 ;
end
if nargin<5 
    dT = 1;
end

N = round(T/dT);
dT = T/N;

PP = zeros(N+1,2);
PP(1,:) = P;

azs = zeros(N,1);
arclens = zeros(N,1);

for ii = 1:N
% [gcdist,p2eaz] = distance(PP(ii,1),PP(ii,2),E(1),E(2),[Re,f]);
[gcdist,p2eaz] = distance(PP(ii,1),PP(ii,2),E(1),E(2));
arclen = V*sind(gcdist)*dT;
[PP(ii+1,1),PP(ii+1,2)] = reckon(PP(ii,1),PP(ii,2),arclen,p2eaz+90);
% convergence issues
% [~,p2eaz2] = distance(PP(ii+1,1),PP(ii+1,2),E(1),E(2));
% [PP(ii+1,1),PP(ii+1,2)] = reckon(PP(ii,1),PP(ii,2),arclen,mean([p2eaz,p2eaz2])+90);

azs(ii) = mod(p2eaz+90,360);
arclens(ii) = arclen;

end
end

