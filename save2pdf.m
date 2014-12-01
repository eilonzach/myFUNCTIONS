function [filepath] = save2pdf(fignumber,filename,odir)
% function [filepath] = save2pdf(fignumber,filename,odir)
% 
% function to save a figure to a good-quality pdf
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

% thisdir = pwd;
% odir = cd(odir);
% cd(thisdir)

if strcmp(odir(end),'/')~=1
    odir = strcat(odir,'/');
end

h = figure(fignumber);

%% Figure size sorting
scrn_wd = 11.28*2.54; % inches to cm
scrn_ht =  7.05*2.54; % inches to cm
pos = get(h,'position');
switch get(h,'Units')
    case 'pixels'
        wd = pos(3)*(scrn_wd/1440); % in cm
        ht = pos(4)*(scrn_ht/800); % in cm
    case 'centimeters'
        wd = pos(3); % in cm
        ht = pos(4); % in cm
end
 
if wd > ht
    orientation = 'landscape';
else
    orientation = 'portrait';
end




set(h,'papertype','usletter');
set(h,'paperorientation',orientation);
set(h,'PaperUnits','centimeters');
set(h,'PaperPosition', [0 0 wd ht]);
print(h,'-painters','-dpdf','-r200',strcat(odir,filename))

filepath = strcat(odir,filename,'.pdf');
end