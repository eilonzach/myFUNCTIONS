function wdir = cd2wdir(fullpath)
%  wdir = cd2wdir(fullpath)
% 
% function to go to the directory in which a certain script is sitting. The
% input is the full path (as a text string) to the  script. If you want to
% use this dynamically to go to whichever directory a certain script you
% are running lives in, type "cd2wdir(mfilename('fullpath'))"

[~,revpath] = strtok(fliplr(fullpath),'/');
wdir = fliplr(revpath);

cd(wdir)

end
