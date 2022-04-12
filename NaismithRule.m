function [time_hr, speed_mph, time_hr_min] = NaismithRule(dist_mi,elevgain_ft)  
% [time_hr, speed_mph, time_hr_min] = NaismithRule(dist_mi,elevgain_ft)  
%  Function to calculate predicted hike duration according to Naismith's
%  rule of 3mi/hr + 30min/1000ft elevation gain
%  INPUTS
%   dist_mi     = distance, in units of miles
%   elevgain_ft = total elevation gain, in units of feet
%  OUTPUTS
%   time_hr     = total time predicted, units of hours
%   speed_mph   = overall speed, units of miles/hr
%   time_hr_min = vector of [hr,min] (rounded up to nearest minute)

    if nargin < 2 || isempty(elevgain_ft)
        elevgain_ft = 0;
    end

    time_hr = dist_mi/3 + 0.5*elevgain_ft/1000;

    speed_mph = dist_mi./time_hr;

    time_hr_min = [floor(time_hr), ceil(60*mod(time_hr,1))];

end