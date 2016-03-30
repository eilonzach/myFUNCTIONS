function [ datf ] = filt_quick( dat,flo,fhi,dt,npoles,filtopt )
%[ datf ] = filt_quick( dat,fmin,fmax,dt )
%   quick function to filter a trace (dat is Nx1 vector) or a set of traces
%   (dat is NxM matrix)

if nargin<4
    dt = 1;
end
if nargin < 5
    npoles = 4;
end
if nargin < 6
    filtopt = 'butter';
end

% calc the filter parms
switch filtopt
    case 'butter'
        % option 1: butter
        [bb,aa]=butter(npoles, [flo, fhi].*dt.*2);
    case 'cheby'
        % option 2: cheby
        [bb,aa]=cheby1(npoles,0.5, [flo, fhi].*dt.*2);
    case 'highbutter'
        % option 3: higher order butter
        [z,p,k]=butter(npoles, [flo, fhi].*dt.*2.);
        [sos,g]=zp2sos(z,p,k); bb=sos; aa=g;
end

% do the filtering
datf = zeros(size(dat));
for ii = size(dat,2)
    datf(:,ii)=filtfilt(bb, aa, dat(:,ii));
end


end

