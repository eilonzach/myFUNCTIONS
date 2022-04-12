function yr = serial2yr(serialtime)
% yr = serial2yr(serialtime)
%  Function to get the year corresponding to a serial time
%  Also works if the input time is a date string

A = datevec(serialtime);
yr = A(:,1);
end