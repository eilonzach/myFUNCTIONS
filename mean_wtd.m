function [mu_w,std_w] = mean_wtd(X,wts)
% [mu_w,std_w] = mean_wtd(X,wts)
% 
% Function to calculate a weighted mean (and, optionally, a standard
% deviation), whereby the mean of the data vector, X, is calculated by
% weighting each datum by the value in the corresponding weights vector,
% wts. 
% Reference: Bevington, P. R., Data Reduction and Error Analysis for
% the Physical Sciences, 336 pp., McGraw-Hill, 1969.

if nargin<2 || isempty(wts) || ~any(wts~=1)
    wts = ones(size(X));
end


totXw = sum(X.*wts);
totw = sum(wts);

% values are complex?
if ~isreal(totXw), totXw = abs(totXw); end
if ~isreal(totw), totw = abs(totw); end

mu_w = totXw./totw;

if nargout > 1 % calculate standard deviation, too
    % effective n of measurements
    n_eff = sum(abs(wts)).^2/sum(abs(wts).^2);
    
    % weightedsum of squared differences
    totdXw = sum( wts.* (X-mu_w).^2 );
    if ~isreal(totdXw), totdXw = abs(totdXw); end
    
    var_w = (totdXw/totw) * (n_eff/(n_eff-1));
    % unbiased std (divide out nu_eff)
    std_w = sqrt(var_w/n_eff);
end

end