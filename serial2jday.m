function jday = serial2jday(serialtime)
% jday = serial2jday(serialtime)
%  Function to get julian day (day of year) corresponding to a serial time
%  Also works if the input time is a date string

A = datevec(serialtime);
jday = doy(A(:,1),A(:,2),A(:,3));

end