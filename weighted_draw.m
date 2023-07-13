function inds = weighted_draw(wts,Ndraw)
% function inds = weighted_draw(wts)
% 
% function to generate a series of random integer indices, "ind", where the
% probability of an index being selected is proportional (linearly) to its
% weight within the vector "wts". I.e. if the weights vector is [1 2 5 2]
% then we expect 10% of the indices to be 1, 20% of the indices to be 2,
% 50% of the indices to be 3, and 20% of the indices to be 4. 

if nargin < 2 || isempty(Ndraw)
    Ndraw = 1; % default is a single weighted random drawn value
end

N = length(wts);
wts = wts(:);
totwt = sum(wts);
pi = [0;cumsum(wts)]./totwt;
inds = ceil(interp1(pi,[0:N],rand(Ndraw,1)));


end
