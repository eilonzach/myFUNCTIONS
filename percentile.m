function val = percentile(vect,pctl,opt)
% val = percentile(vect,pctl,opt)
%
% Function to get the pctl percentile value of a series of sorted values,
% arranged from smallest to largest. 
% opt can take the values 'interp' (the default) or 'roundnear' to
% choose whether to interpolate the sorted values to get the nth
% percentile, or whether to simply take the sorted value with index closest
% to the chosen percentile*the data length

if nargin < 3
    opt = 'interp';
end

N = length(vect);
vects = sort(vect);
if strcmp(opt,'interp')
    ivect = linspace(0,100,N);
    val = interp1(ivect,vects,pctl);
elseif strcmp(opt,'roundnear')
    val = vects(round(N.*pctl/100));
end

end