function [dof,samp_per_fp]=scdofcalc(xx,option)
% [dof,tim]=scdofcalc(xx,option)
% scdofcalc.m  calculate degrees of freedom a la Silver and Chan for a time
% series xx or from a noise series XY filename
% optional arguments: 
% [] = xx is the time series. 
% 'filenam' = xx is the name of the file with the time series

if nargin==2
    if strcmp(option,'filenam')==1
        infile=xx;
        timeseries=load(infile);
    end
else
    timeseries=xx;
end

f=fft(timeseries(:,1)); %NB was f=fft(timeseries(:,2));
f2=f.*conj(f);
nyq=floor(length(f2)/2);
f2 = f2(1:nyq);
f4 = f2 .* f2;
ee = sum(f2) - 0.5*(f2(1) + f2(nyq));
e4 = sum(f4) - 0.5*(f4(1) + f4(nyq));

dof =  2 * (2*ee*ee/e4 - 1);		% SC equation A12
samp_per_fp=length(timeseries)/dof; 




ifplot = false; if ifplot; warning('Setting ifplot == 1'); end 
if ifplot; 
    fprintf('\nDegrees of freedom: %1.2f. nsamp = %1.0f. samp per free parameter??? = %3.0f. \n',...
        dof, length(timeseries), samp_per_fp); 

    figure(1111); clf; hold on; 
    tiledlayout(2,1,'TileSpacing', 'compact'); 
    nexttile(); hold on; box on; grid on; 
    title(sprintf('Time series. N = %1.0f',length(timeseries))); 
    xlabel('Index'); 
    plot(timeseries); 
    nexttile(); hold on; box on; grid on; 
    title(sprintf('Autocorrelation. DOF = %1.4f', dof)); 
    xlabel('Index'); 
    [corr, lags] = xcorr(timeseries); 
    corr = corr ./ max(corr); 
    plot(lags, corr); 
end




end