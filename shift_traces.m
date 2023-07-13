function odat = shift_traces(idat,tshift,dt,opt)
% odat = shift_traces(idat,tshift,dt)
% 
% Function to shift traces by some time.
% 
% Default is for second argument to be time in seconds, third argument to
% be data dt values. If no third argument, then assumes tshift is actually
% an integer number of samples to shift.
% 
% POSITIVE values of tshift correspond to DELAYS to be imposed. I.e. a
% positive tshift will move the features in the time series later in time.
% 
% If idat contains multiple columns, tshift should/can be a vector of
% different shifts for each column. 
% 
% "opt" allows user to specify whether time shift should be "rounded"
% (default) whereby shift is by rounded integer number of samples, or
% "precise", whereby shift is done by interpolation. 
% 
% Note, shifting will add samples to one end of the time series, and
% subtract samples from the other end. The added samples will by default be
% zeros. It is recommended that idat be tapered before shifting. 


if nargin < 3 || isempty(dt)
    dt = 1; 
end
if nargin < 4 || isempty(opt)
    opt = 'rounded'; 
end

Nsamp = size(idat,1);
Ntrace = size(idat,2);
odat = zeros(Nsamp,Ntrace);

for ii = 1:Ntrace
    if strcmp(opt,'rounded')
        lagshift = round(tshift(ii)/dt);
        if lagshift >=0
            odat(lagshift+1:end,ii) = idat(1:end-lagshift,ii);
        else
            odat(1:end+lagshift,ii) = idat(1-lagshift:end,ii);
        end
    elseif strcmp(opt,'precise')
        tt0 = dt*[0:(Nsamp-1)]';
        tt1 = dt*[0:(Nsamp-1)]' - tshift; 
        odat(:,ii) = interp1(tt0,idat(:,ii),tt1,"linear",0);
    else
        error('opt not correctly specified')
    end
end

