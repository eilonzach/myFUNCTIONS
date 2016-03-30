function [MM,DD] = calday(YYYY,JJJ,mode)
%[MM,DD] = calday(YYYY,JJJ,mode)
% converts day of year into calendar date (month & day)
% if mode is 's', will output two 2-character strings
% if mode is 'i', will output two numbers
%
% Written by Zach Eilon, 2011


if nargin==1
    JJJ=YYYY;
    YYYY=2010;
end

A=ones(31,12);
dpm=eomday(YYYY,[1:12]);
for k=1:length(dpm)
    A(dpm(k)+1:end,k)=0;
end
A(A>0)=[1:sum(dpm)]';
[DD,MM]=find(A==JJJ);

if nargin==3
    if strcmp(mode,'s')==1
        if MM >= 10
            MM=sprintf('%s',num2str(MM));
        elseif MM < 10 
            MM=sprintf('0%s',num2str(MM));
        end
        if DD >= 10
            DD=sprintf('%s',num2str(DD));
        elseif DD < 10 
            DD=sprintf('0%s',num2str(DD));
        end 
    elseif strcmp(mode,'i')==0
        error('\tMode "%s" not recognised',mode);
        return
    end
end
end

