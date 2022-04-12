function y_ = rm_intercept(y,x,wts)
% y_ = rm_intercept(y,x,wts)
% 
% function to remove y intercept from a linear trend - fit by least squares
% if no x values, assume x = 0:length(y)-1; Can weight data by wts.

y = y(:); % ensure column. 
N = length(y);

if nargin < 2 || isempty(x)
    x = [0:N-1]';
end
if nargin < 3 || isempty(wts)
    wts = ones(N,1);
end



% make G for linear fit - first model parm is intercept;
G = [ones(N,1),x];
iCd = diag(wts);

m_est = (G'*iCd*G)\G'*iCd*y;

b = m_est(1);

y_ = y-b;
    


end