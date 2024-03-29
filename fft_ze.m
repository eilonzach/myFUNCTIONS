function [ XX,ff,XX_0,ff_0 ] = fft_ze( xx,tt )
%[ XX,ff,XX_0,ff_0 ] = fft_ze( xx,tt )
%   Function to do fast fourier transform of a data series, xx, to also get
%   the frequencies. Second two outputs make the spectrum symmetric about
%   zero frequency - better for plotting
% 
%  INPUTS:
%    xx     time series vector of Nx1 points (or NxM) if M traces
%    tt     time vector (or, if scalar then just dt)
% 
%  OUTPUTS:
%    XX     fast fourier transform of xx (a complex, Nx1 vector)
%    ff     frequencies corresponding to XX
%    XX_0   XX but made symmetric about freq=0
%    ff_0   ff but made symmetric about freq=0

if nargin < 2
    tt = 1;
end

if any(size(xx)==1)
    xx = xx(:); % columnise if vector. If matrix, assume already in columns
end

%% get dt
if isscalar(tt)
    dt = tt;
else
    dt = mean(diff(tt));
end

%% get important vals
N = size(xx,1);
T = dt*N;


%% make freq. axis
if mod(N,2)
     ff = [0:(N-1)/2,-(N-1)/2:-1]*(1/T);
else
     ff = [0:N/2,-N/2+1:-1]*(1/T);
end
ff = ff(:);

%% add in zero-padding for stability

%% do fft
XX = fft(xx,[],1);

%% make 0-symmetric versions
ind = [ceil(N/2)+1:N,1:ceil(N/2)];
XX_0 = XX(ind,:);
ff_0 = ff(ind);


end

