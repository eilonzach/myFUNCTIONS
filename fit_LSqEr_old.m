function [ m,b,mstd,bstd,m_upper,m_lower,misfit_best,R2,Covmat,wbmisfitratio ] = fit_LSqEr_old( X,Y,b0opt,normopt,mprior,bprior,xwt,ywt )
%[ m,b,mstd,bstd,m_upper,m_lower,misfit_best ,R2,Covm,wbmisfitratio] = fit_LSqEr_old( X,Y,b0opt=1,normopt=2,mprior,bprior,xwt,ywt )
%   function to fit a straight line to X,Y data where there is error in
%   both data by grid searching to find the line that gives the minimum
%   squared misfit in both dimensions with a line equation y = m*x + b. 
% 
% upper-lower are the 1-sigma bounds, as dictated by F-test
%   
% The final two outputs are the covariance matrix (Covmat) and the ratio
% between the misfits calculated for the preferred model, and the model
% with a gradient 90deg rotated from the preferred (if a cloud of points,
% not v. different. If a linear array, hugely different)

if nargin<3 || isempty(b0opt)
    b0opt = 1; %option to pin b to zero
end
if nargin<4 || isempty(normopt)
    normopt = 2;
end
if nargin < 5 || isempty(mprior)
    % find something about the bounds of the data
    dx = max(X)-min(X);
    dy = max(Y)-min(Y);
    m0 = dy/dx;
    mprior = -1000:1000;
end
if nargin < 6 || isempty(bprior)
    if b0opt
        bprior = 0;
    else
        bprior = -1000:1000;
    end
end
if nargin < 7
    xwt = ones(size(X));
end
if nargin < 8
    ywt = ones(size(X));
end

M = length(mprior);
B = length(bprior);

misfit = zeros(M,B);

for im = 1:length(mprior)
    for ib = 1:length(bprior)
        dy = Y - ( mprior(im)*X + bprior(ib) );
        dx = X - ( (Y-bprior(ib))./mprior(im) );
        dy = ywt.*dy;
        dx = xwt.*dx;
        misfit(im,ib) = norm(dy,normopt) + norm(dx,normopt);
    end
end

[misfit_best,ib_best,im_best] = mingrid(misfit);

m = mprior(im_best);
b = bprior(ib_best);

%% F tests to get confidence bounds for m
nar = length(X);
Fval_m = zeros(M,1);
for ii = 1:M
    Fval_m(ii) = (misfit(im_best,ib_best)\misfit(ii,ib_best));
end
Fcrit = finv(0.67,nar-1,nar-1); % for 1 std
m_upper = mprior(find(Fval_m<Fcrit,1,'last'));
m_lower = mprior(find(Fval_m<Fcrit,1,'first'));

%% F tests to get confidence bounds for b
nar = length(X);
Fval_b = zeros(B,1);
for ii = 1:B
    Fval_b(ii) = (misfit(im_best,ib_best)\misfit(im_best,ii));
end
Fcrit = finv(0.67,nar-1,nar-1); % for 1 std
b_upper = bprior(find(Fval_b<Fcrit,1,'last'));
b_lower = bprior(find(Fval_b<Fcrit,1,'first'));

mstd = (m_upper-m_lower)/2;
bstd = (b_upper-b_lower)/2;

%% R-squared
% R2 = 1 - SSE/SST (unadjusted)

SSE = sum( ywt.*(Y-(m*X+b)).^2 ) + sum( xwt.*(X-(Y-b)/m).^2 );
SST = sum( ywt.*Y.^2 ) + sum( xwt.*X.^2 );

R2 = 1 - SSE/SST;

%% covariance
if nargout >8
    Covmat = cov(X,Y);
end

%% misfit ratio
% have to demean the data in order to do this 90deg flip fit
dy_ = demean(Y,'constant') - ( -demean(X,'constant')./m );
dx_ = demean(X,'constant') - ( -demean(Y,'constant')*m );
dy_ = ywt.*dy_;
dx_ = xwt.*dx_;
misfit_worst = norm(dy_,normopt) + norm(dx_,normopt);
wbmisfitratio = misfit_worst./misfit_best;

%% p-value
% need to define the null hypothesis. Null state for sampling two random,
% and unconnected variables, is a joint pdf with semi-minor and -major axes
% parallel to the two coordinate axes. Or, an isotropic cloud of points. 


end

