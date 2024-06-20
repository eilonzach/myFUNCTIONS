function [coeffs,lambdas] = fourier_decompose(yy,xx,minmaxlambda,ifplot)
% decompose a series into its sinusoidal components
% output is Nk*2 + 1 coefficients, 
% first Nk are cosine terms
% next Nk are sine terms
% minmaxlambda are limiting wavelengths to consider (same units as xx) and
%   if only one is given, then assume max.
% last coeff is constant term (mean)

% plan for nan
nnany = find(~isnan(yy));

% find wavenumbers, build frequency grids
maxk = 2*pi./mean(diff(xx))/2;
mink = 2*pi./(xx(end)-xx(1));
kk = [mink:mink:maxk];
lambdas = 2*pi./kk;

if ~isempty(minmaxlambda)
    if length(minmaxlambda)==1
        oklambda = lambdas<=minmaxlambda;
    elseif length(minmaxlambda)==2
        oklambda = lambdas<=minmaxlambda(2) & lambdas >= minmaxlambda(1);
    end
    kk = kk(oklambda);
    lambdas = lambdas(oklambda);
end


Nk = length(kk);
Nx = length(xx);

xx_orig = xx;
xx = xx(nnany);


G = make_sincos_G(xx,kk);
if ifplot
G_orig = make_sincos_G(xx_orig,kk);
end
% LSQR
coeffs = (G'*G)\G'*yy(nnany);

% RMS error
rms(G*coeffs - yy(nnany));

if ifplot
    figure(154), clf
    plot(xx,yy(nnany),'o-k','linewidth',3)
    hold on
    % cos terms
    plot(xx_orig,G_orig(:,1:Nk).*(ones(Nx,1)*coeffs(1:Nk)') + coeffs(end),...
        'c','linewidth',0.5);
    % sin terms
    plot(xx_orig,G_orig(:,Nk+[1:Nk]).*(ones(Nx,1)*coeffs(Nk+[1:Nk])') + coeffs(end),...
        'b','linewidth',0.5);
    % final fit
    plot(xx_orig,G_orig*coeffs,'r--','linewidth',2.5)
end

end

function G = make_sincos_G(xx,kk)
    Nk = length(kk);
    Nx = length(xx);
    % build G matrix for LSQR fit
    G = zeros(Nx,Nk*2+1);
    G(:,end) = 1; % constant term at end
    for ik = 1:Nk
        G(:,ik) = cos(kk(ik)*xx);
        G(:,ik+Nk) = sin(kk(ik)*xx);
    end
end