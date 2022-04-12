function lolaz = geo_UTM2lalo(enz,lon_zone,hemi)
% ofile = geo_UTM2lalo(ifile,zone,hemi)
%
% Function to convert UTM* coordinates into geographic coordinates
% 
% INPUTS:
%   1) an .xyz file or matrix with two(three) columns:  
%      Easting, Northing, (and Depth), with E and N  in UTM* coordinates
%   2) the longitudinal zone, either as a string (with optional
%      Latitude band, eg. '32H') or as a decimal value equal to the centre
%      longitude of the data (from which  function will compute the zone)
%   3) Hemisphere, 'N'|'North'|+1 or 'S'|'South'|-1. Default is N.
% 
% OUTPUTS:
%     a .lola(z) file or a [lon lat (z)] matrix, depending on input type
% 
% 
%  * Universtal Transverse Mercator system
%  https://en.wikipedia.org/wiki/Universal_Transverse_Mercator_coordinate_system
% 
%  Z. Eilon 3/2021

if ischar(enz) % input is a file
    otype = 'file';
    ifile = enz;
    enz = load(ifile);
else
    otype = 'matrix';
end

ncol = size(enz,2);

if ncol < 1 || ncol > 4
    error('This file does not appear to have 2 or 3 columns')
end

E = enz(:,1);
N = enz(:,2);

if ischar(lon_zone) % lon zone is a code
    % convert from string with numbers in code
    try % first assume code is two character number + 1 char letter (which we ignore)
        lonzone = eval(lon_zone(1:2));
    catch % ok so maybe code is 1 character number + 1 char letter (which we ignore)
        lonzone = eval(lon_zone(1));
    end
    lon0 = lonzone*6 - 183;
else % lon zone is a longitude value of the centre of the array
    % need to round to centre of bin
    lon0 = (fix(mod(lon_zone+180,360)/6)+1)*6 - 183;
end

ns = 1;
if any(strcmp({'S','s','South','south'},hemi)) || hemi==-1
    ns = -1;
end

%% intermediate values (lengths in *km*)
% put input data into km
if any(abs(N)>10000) % in m
    N = N/1e3;
end
if any(abs(E)>20000) % in m
    E = E/1e3;
end

if ns == 1
    N0 = 0;
elseif ns == -1
    N0 = 10000;
end
E0 = 500;
k0 = 0.9996;
f = 1/298.257223563; % flattening
a = 6378.137; % equatorial radius

n = f/(2-f);
n2 = n.^2;
n3 = n.^3;
n4 = n.^4;

A = (a/(1+n))*(1 + n2/4 + n4/64); %not clear what this series is supposed to be

alp(1,1) = n/2 - 2*n2/3 + 5*n3/16; 
alp(1,2) = 13*n2/48 - 3*n3/5;
alp(1,3) = 61*n3/240;
bet(1,1) = n/2 - 2*n2/3 + 37*n3/96;
bet(1,2) = n2/48 + n3/15; 
bet(1,3) = 17*n3/480;
del(1,1) = 2*n - 2*n2/3 - 2*n3;
del(1,2) = 7*n2/3 - 8*n3/5;
del(1,3) = 56*n3/15;

xi = (N - N0)/(k0*A);
eta = (E - E0)/(k0*A);

j = [1,2,3];
ndat = length(N);
xi_ = xi - sum((ones(ndat,1)*bet).*sin(2*xi*j).*cosh(2*eta*j),2);
eta_ = eta - sum((ones(ndat,1)*bet).*cos(2*xi*j).*sinh(2*eta*j),2);
% sig_ = 1 - 2*sum((ones(ndat,1)*(j.*bet)).*cos(2*xi*j).*cosh(2*eta*j),2);
% tau_ = 2*sum((ones(ndat,1)*(j.*bet)).*sin(2*xi*j).*sinh(2*eta*j),2);

chi = asind(sin(xi_)./cosh(eta_)); % degrees now

% now actual lat/lon
lat = chi + sum((ones(ndat,1)*del).*sind(2*chi*j),2); % degrees now
lon = lon0 + atand(sinh(eta_)./cos(xi_)); % degrees

enz(:,1) = lon;
enz(:,2) = lat;

switch otype
    case 'file'
        ofile = [ifile,'.lolaz'];
        ofile = regexprep(ofile,'.xyz','');
        fid = fopen(ofile,'w');
        for ii = 1:length(lon)
            fprintf(fid,'%10.5f %10.5f %f\n',enz(ii,:));
        end
        fclose(fid);
        lolaz = ofile;
    case 'matrix'
        lolaz = enz;
end


end