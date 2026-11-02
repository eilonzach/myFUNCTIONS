function dt = serial2datetime(sr)
% dt = serial2datetime(sr)
%
% function to turn a serial time value (probably ~7e5, in days since Jan 0)
% into a matlab datetime value

dt = datetime(datevec(sr));

end
