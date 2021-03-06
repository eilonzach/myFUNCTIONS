function [filepath] = save2png_basic(fignumber,filename,odir)
% function [filepath] = save2png_basic(fignumber,filename,odir)
% 
% function to save a figure to a basic-quality png
%
% INPUTS:
% fignumber - (required) number of figure to save
% filename - (optional, default is figN) name of saved file
% odir - (optional, default is current dir) dir to save file in

if nargin < 3
    odir = pwd;
end
if nargin < 2
    filename = sprintf('fig%u',fignumber);
end

if ~strcmp(filename(end-3:end),'.png')
    filename = [filename '.png'];
end

if strcmp(odir(end),'/')~=1
    odir = strcat(odir,'/');
end

h = figure(fignumber);
print(h,'-dpng',strcat(odir,filename))

filepath = strcat(odir,filename);
end