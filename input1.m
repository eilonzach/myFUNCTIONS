function [ out ] = input1(prompt,numopt)
% function to get one character of user input without pressing enter
%
% 'prompt' is the prompt that appears at the command line
% 'stropt' is an option to specify that the input is a number, not a string

if nargin<2
    numopt=0;
end

fprintf(prompt)

w = waitforbuttonpress;
if w 
    out = get(gcf,'CurrentCharacter');
end

if numopt==1
    out = str2double(out);
end

end

