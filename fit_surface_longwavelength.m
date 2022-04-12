function [dat,m_sin,m_lin,results] = fit_surface_longwavelength(dat,x_wvls,y_wvls,rmlin,cliplevel,ifplot,ndownsamp,datastr)
% [dat,m_sin,m_lin,results] = fit_surface_longwavelength(dat,x_wvls,y_wvls,rmlin,cliplevel,ifplot,ndownsamp)
%
% function to fit surface with long-wavelength 2D structure. Assumes you
% already have structure dat, which contains fields 'x','y','z', in
% cartesian coordinate system, x/y distances in km.

% default parms - set if not specified
if nargin < 2 || (isempty(x_wvls) && isempty(y_wvls))
    x_wvls = [200:200:400,500,1200];
    y_wvls = [200:200:400,500,1200];
end
if nargin < 4 || isempty(rmlin)
    rmlin = true;
end
if nargin < 5 || isempty(cliplevel)
    cliplevel = 10; % percent on either end of distribution to clip off
end
if nargin < 6 || isempty(ifplot) 
    ifplot = 0;
end
if nargin < 7 || isempty(ndownsamp) 
    ndownsamp = 40;
end
if nargin < 8 || isempty(datastr) 
    datastr = 'data';
end

clip_or_nan = 'nan'; % if 'clip' will set as max/min values. If 'nan' will delete points
if ~isfield(dat,'Npts')
    dat.Npts = length(dat.z);
end

% check column vectors
dat.x = dat.x(:);
dat.y = dat.y(:);
dat.z = dat.z(:);

symsize = 40;

%% clip extremes of data
sortz = sort(dat.z(~isnan(dat.z)));
minz = sortz(max([round((cliplevel/100)*length(sortz)),1]));
maxz = sortz(round((1 - cliplevel/100)*length(sortz)));
meanz = mean(dat.z(dat.z <= maxz & dat.z >= minz)); % get mean value, excluding outliers

% make clipped, de-meaned depths
dat.z_cd = dat.z;
if strcmp(clip_or_nan,'clip')
    dat.z_cd(dat.z_cd > maxz) = maxz; % N.B. choice is to clip values outside cliplevel, not NaN them
    dat.z_cd(dat.z_cd < minz) = minz;
elseif strcmp(clip_or_nan,'nan')
    dat.z_cd(dat.z_cd > maxz) = nan; % N.B. choice is to clip values outside cliplevel, not NaN them
    dat.z_cd(dat.z_cd < minz) = nan;
end
dat.z_cd = dat.z_cd - meanz;

% get rid of any extreme x,y points (outside main array)
% - no don't do this right now -

%% get rid of any nan data 
nnan = ~isnan(dat.z_cd);
% save originals
dat.x_orig = dat.x;
dat.y_orig = dat.y;
dat.z_orig = dat.z;
% delete nan points
dat.x = dat.x(nnan);
dat.y = dat.y(nnan);
dat.z = dat.z(nnan);
dat.z_cd = dat.z_cd(nnan);
dat.Npts = length(dat.z);

%% downsampler to make plotting not take so long...
downsamp = [1:ndownsamp:dat.Npts];

% while testing
% figure(9); plot(sort(dat.z_cd))
% testing

% while testing
% figure(11), clf
% subplot(121)
% scatter(dat.lon(downsamp),dat.lat(downsamp),20,dat.z_cd(downsamp))
% 
% subplot(122)
% scatter(dat.x(downsamp),dat.y(downsamp),20,dat.z_cd(downsamp))
% colorbar

%% solve for best fitting linear trend
% m = coefficients of [x,y,xy,-xy]
G_lin = [ones(size(dat.x)),dat.x,dat.y,dat.x.*dat.y,-dat.x.*dat.y];
% G_lin = [ones(size(dat.x)),dat.x,dat.y];

[m_lin,flag,relres,iter,resvec] = lsqr( G_lin, dat.z_cd, 1e-5, 500 );
if iter==500, fprintf('Warning LSQR not converging\n'); end

% compute predicted depths for linear trends
z_lin = G_lin*m_lin;
% compute residual from linear fit
dat.z_cd_nl = dat.z_cd - z_lin;

% consider misfit
% figure(13)
% plot(dat.z_cd,z_lin,'o')

if ifplot
    figure(ifplot), clf, set(gcf,'pos',[40 142 2093 1121]);
    subplot(231)
    scatter(dat.x(downsamp),dat.y(downsamp),symsize,dat.z_cd(downsamp),'filled')
    xlabel('SE ==>   (x)','fontsize',18,'fontweight','bold')
    ylabel('NE ==>   (y)','fontsize',18,'fontweight','bold')
    title(['De-meaned, clipped ',datastr],'fontsize',20,'fontweight','bold')
    caxis([minz maxz])

    subplot(232)
    scatter(dat.x(downsamp),dat.y(downsamp),10,z_lin(downsamp),'filled')
    xlabel('SE ==>   (x)','fontsize',18,'fontweight','bold')
    ylabel('NE ==>   (y)','fontsize',18,'fontweight','bold')
    title(['Linear predicted ',datastr],'fontsize',20,'fontweight','bold')
    caxis([minz maxz])

    subplot(233)
    scatter(dat.x(downsamp),dat.y(downsamp),symsize,dat.z_cd_nl(downsamp),'filled')
    xlabel('SE ==>   (x)','fontsize',18,'fontweight','bold')
    ylabel('NE ==>   (y)','fontsize',18,'fontweight','bold')
    title(['Non-linear-fit ',datastr],'fontsize',20,'fontweight','bold')
    caxis([minz maxz])
    colorbar
end

%% solve for best fitting sinusoidal terms
% set up results table
results = table('size',[length(x_wvls)+length(y_wvls),3],'VariableNames',{'Amplitude','Offset','Lambda'},'VariableTypes',{'double','double','double'});
G_sin = [];
irow = 0;
for iwx = 1:length(x_wvls)
    G_sin = [G_sin, sin(2*pi*dat.x./x_wvls(iwx)), cos(2*pi*dat.x./x_wvls(iwx))];
    irow = irow+1;
    results.Properties.RowNames{irow} = sprintf('L_x=%.0fkm',x_wvls(iwx));
    results.Lambda(irow) = x_wvls(iwx);
end
for iwy = 1:length(y_wvls)
    G_sin = [G_sin, sin(2*pi*dat.y./y_wvls(iwy)), cos(2*pi*dat.y./y_wvls(iwy))];
    irow = irow+1;
    results.Properties.RowNames{irow} = sprintf('L_y=%.0fkm',y_wvls(iwy));
    results.Lambda(irow) = y_wvls(iwy);
end
M_sin = size(G_sin,2);

if rmlin
    z_sin_toinv = dat.z_cd_nl; % dat.z_cd or dat.z_cd_nl
else
    z_sin_toinv = dat.z_cd; % dat.z_cd or dat.z_cd_nl
end

[m_sin,flag,relres,iter,resvec] = lsqr( G_sin, z_sin_toinv, 1e-5, 500 );
% parse results into coefficients
for iw = 1:length(m_sin)/2
    s = m_sin(2*iw - 1);
    c = m_sin(2*iw);
    results{iw,1} = round(sqrt(s.^2 + c.^2),3);
    results{iw,2} = round(atan2(c,s)*results.Lambda(iw)/2/pi,1);
end
results

if ifplot
    figure(9), set(gcf,'pos',[1979 890 560 420])
    plot(results.Lambda,results.Amplitude,'-ok','linewidth',2)
    xlabel('Lambda (km)'); ylabel('Amplitude');
end

% compute predicted depths 
if rmlin
    z_pred = G_sin*m_sin + G_lin*m_lin;
else
    z_pred = G_sin*m_sin;
end
dat.z_pred = z_pred;

% corrected/residual bathymetry
dat.z_unfit = dat.z_cd - z_pred;

if ifplot
    figure(ifplot)
    subplot(235)
    scatter(dat.x(downsamp),dat.y(downsamp),10,z_pred(downsamp),'filled')
    xlabel('SE ==>   (x)','fontsize',18,'fontweight','bold')
    ylabel('NE ==>   (y)','fontsize',18,'fontweight','bold')
    title(['Sinusoidal predicted ',datastr],'fontsize',20,'fontweight','bold')
    caxis([minz maxz])
    colorbar

    subplot(236)
    scatter(dat.x(downsamp),dat.y(downsamp),symsize,dat.z_unfit(downsamp),'filled')
    xlabel('SE ==>   (x)','fontsize',18,'fontweight','bold')
    ylabel('NE ==>   (y)','fontsize',18,'fontweight','bold')
    title(['Residual (un-fit) ',datastr],'fontsize',20,'fontweight','bold')
    caxis([minz maxz])
    colorbar
end

%% norm reductions
fprintf('Norm of original %s data: \t\t%10.1f km\n',datastr,norm(dat.z))
fprintf('Norm of clipped, demeaned %s data: \t%10.1f km\n',datastr,norm(dat.z_cd))
fprintf('Norm of non-linear-fit %s data: \t\t%10.1f km\n',datastr,norm(dat.z_cd_nl))
fprintf('Norm of non-sin-fit %s data: \t\t%10.1f km\n',datastr,norm(dat.z_unfit))

%% variance reductions (compared to clipped, demeaned)
addpath('/Users/zeilon/Dropbox/MATLAB/BWTOMOG_atten_Vp_Vs_joint/matguts/')
fprintf('VR for linear-fit %s data: %4.1f %%\n',datastr,variance_reduction(dat.z_cd,z_lin))
fprintf('VR for sin-fit %s data:    %4.1f %%\n',datastr,variance_reduction(dat.z_cd,z_pred))

% plot predictions
% figure(14); hold on
% plot(dat.z_cd,z_lin,'o','Markersize',2)
% plot(dat.z_cd,z_pred,'o','Markersize',2)

%% Y-components of bathymetry (best fitting 2D)
isy = 2*length(x_wvls) + [1:2*length(y_wvls)];
yy = [min(dat.y):1:max(dat.y)]';
Gy_sin = [];
for iwy = 1:length(y_wvls)
    Gy_sin = [Gy_sin, sin(2*pi*yy./y_wvls(iwy)), cos(2*pi*yy./y_wvls(iwy))];
end
% account for x-offset at x=0 (cosine terms)
zx = sum(m_sin(2:2:2*length(x_wvls)));
if rmlin
    dat.zy = Gy_sin*m_sin(isy) + yy.*m_lin(3) + meanz + zx; % dat.z_cd or dat.z_cd_nl
else
    dat.zy = Gy_sin*m_sin(isy) + meanz + zx; % dat.z_cd or dat.z_cd_nl
end

if ifplot
    subplot(235), hold on
    plot(zeros(size(yy)),yy,'-k','linewidth',2.5)


    subplot(234), cla, hold on
    plot(dat.y(abs(dat.x) <= 10),dat.z(abs(dat.x) <= 10),'ok','markersize',4)
    plot(yy,dat.zy,'r','linewidth',3)
    legend('Z_{obs} within |x|<10','Z(y) predicted','location','southeast')
    set(gca,'xlim',[min(dat.y),max(dat.y)],'fontsize',16,'linewidth',1.5,'layer','top')
    xlabel('MM-parallel (Y) distance (km)','fontsize',18,'fontweight','bold')
    ylabel('Z (m) predicted','fontsize',18,'fontweight','bold')
    title('y-fit only (mean added back in)','fontsize',20,'fontweight','bold')
end

%% outputs
% dat