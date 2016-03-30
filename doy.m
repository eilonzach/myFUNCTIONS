function [JJJ] = doy(YYYY,MM,DD,mode)
%[JJJ] = doy(YYYY,MM,DD,mode)
% converts calendar date (day + month) into day of year
% if mode is 's', will output a 3-character string
% if mode is 'i', will output a number
% 
% Written by Zach Eilon, 2011

if nargin==3
    JJJ=datenum(YYYY,MM,DD)-datenum(YYYY,0,0);
elseif nargin==4
    if strcmp(mode,'s')==1
        JJJ=datenum(YYYY,MM,DD)-datenum(YYYY,0,0);
        if JJJ >= 100
            JJJ=sprintf('%s',num2str(JJJ));
        elseif JJJ >= 10 && JJJ < 100
            JJJ=sprintf('0%s',num2str(JJJ));
        elseif JJJ< 10
            JJJ=sprintf('00%s',num2str(JJJ));
        end
    elseif strcmp(mode,'i')==1
        JJJ=datenum(YYYY,MM,DD)-datenum(YYYY,0,0);
    end
end
end

