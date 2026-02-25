% plot of wavetrain with measured and anticipated parallel magnetic field
% load data from art2.load_data

%% Collect all the parameters
% 1. define time intervals
% 2. make matching/alternatively find these from sometime before
% 3. make fit of to model electric field/potential
% 4. calculate the magnetic field from these parameters
% 5. collect the magnetic field in one timer series

time = {toepoch([2007 08 31 10 17 45.60; 2007 08 31 10 17 45.80])',...
        toepoch([2007 08 31 10 17 45.90; 2007 08 31 10 17 46.08])',...
        toepoch([2007 08 31 10 17 46.08; 2007 08 31 10 17 46.20])',...
        toepoch([2007 08 31 10 17 46.20; 2007 08 31 10 17 46.38])',...
        toepoch([2007 08 31 10 17 46.35; 2007 08 31 10 17 46.46])',...
        toepoch([2007 08 31 10 17 46.46; 2007 08 31 10 17 46.60])'};
doMatch = 0;
if doMatch    
    for kk = 6:numel(time);
        tint = time{kk};
        art2.domatch;
        cn.print(['matching_' num2str(kk)])
    end
end

% the collected parameters, from plots in Days/2014-04-24
phi = {170,120,190,300,170,200};
length = [13,11,15,25,23,10]; % peak to peak
width = length*1.5; % make them spherical, or pancakeish
velocity = {373,419,770,916,938,559};


%% Calculate magnetic field from above parameters       
for kk = 1:numel(width)
    B0 = 25; % nT
    n = 0.04; % cc
    lr = width(kk)/2;
    lz = length(kk)/2;
    nz = 50;
    phi0 = phi{kk};
    nlz = 10;
    art2.calculateB_rz; 
    BZ{kk} = squeeze(Bz)*1e9; % nT
    distance{kk} = zz;
end

%% Collect fields together and make time series
% make length into time
BB = [];
for kk = 1:numel(width)
    timeseries{kk} = time{kk}(1)+(distance{kk}-distance{kk}(1))/velocity{kk};
    BB = [BB; tocolumn(timeseries{kk}) tocolumn(BZ{kk})];
end
BB3 = BB;
BB4 = BB;
%% fine tune magnetic field placements
% do this after having done the figure to get the times of the center of 
% the electron hole
tzoom{1} = toepoch([2007 08 31 10 17 45.3;2007 08 31 10 17 46.9])'; % ind=5, quality=2
%BB3 = [tocolumn(linspace(tzoom{1}(1),tzoom{1}(2),10000)) zeros(10000,1)];
%BB4 = [tocolumn(linspace(tzoom{1}(1),tzoom{1}(2),10000)) zeros(10000,1)];
doAdd = 0;
doPad = 1;
doLoadTime = 1;
doPicktime = 0;

if doLoadTime
    load /Users/Cecilia/Data/BM/20070831/tcenter.mat
elseif doPicktime
    [ts3,~] = get_time(6);
    [ts4,~] = get_time(6);
end
if doAdd; 
    BB3 = []; BB4 = [];
elseif doPad; 
    BB3 = [tocolumn(linspace(tzoom{1}(1),tzoom{1}(2),10000)) zeros(10000,1)];
    BB4 = [tocolumn(linspace(tzoom{1}(1),tzoom{1}(2),10000)) zeros(10000,1)];
end
for kk = 1:numel(width)
    timeseries3{kk} = ts3(kk)+distance{kk}/velocity{kk};
    timeseries4{kk} = ts4(kk)+distance{kk}/velocity{kk};
    newB3 = [tocolumn(timeseries3{kk}) tocolumn(BZ{kk})];        
    newB4 = [tocolumn(timeseries4{kk}) tocolumn(BZ{kk})];        
    if doAdd        
        BB3 = [BB3; newB3];
        BB4 = [BB4; newB4];
    elseif doPad
        newB3 = cn.zeropad(newB3,tzoom{1});
        newB4 = cn.zeropad(newB4,tzoom{1});
        BB3 = irf_add(1,newB3,1,BB3);
        BB4 = irf_add(1,newB4,1,BB4);
    end
    %BB3 = irf_add(1,BB3,1,newB);
    %BB3 = [BB3; tocolumn(timeseries{kk}) tocolumn(BZ{kk})];
        
end
if 0
for kk = 1:numel(width)
    timeseries{kk} = ts4(kk)+distance{kk}/velocity{kk};
    newB = [tocolumn(timeseries{kk}) tocolumn(BZ{kk})];
    %BB4 = irf_add(1,BB4,1,newB);
    %[BB4; tocolumn(timeseries{kk}) tocolumn(BZ{kk})];
    %BB4 = [BB4; irf_tlim(newB,time{kk})];
    BB4 = [BB4; newB];
end
end
    
%% Make plot !!!
% This script produces a figure with either one or two panels with
% wavetrains, that serve to illustrate the high probability that we observe
% the same holes on both sc.
printTitle = 0;
% fig4: a wave train
cd /Users/Cecilia/Data/BM/20070831
tzoom{1} = toepoch([2007 08 31 10 17 45.3;2007 08 31 10 17 46.9])'; % ind=5, quality=2
tzoom{2} = toepoch([2007 08 31 10 17 51.2;2007 08 31 10 17 52.4])';
fmin = 2;
fmax = 0;
fs = cn.f(staBpar3);
order = 5;
BB3 = irf_resamp(BB3,staBpar3);
BB4 = irf_resamp(BB4,staBpar4);
% Make figure
set(0,'defaultLineLineWidth', 1);
fig = figure(24);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');

nPanels = 3;
%set(gcf,'position',[10   684   591   295]);
set(gcf,'position',[10 548 591 408]);
h = irf_plot(nPanels);
isub = 1;
colororder = [0.2 0.5 0.2;0 0 1;0.2 0.5 0.2;0 0 1];
colororderLight = [0.8 1 0.8;0.8 0.8 1;0.2 0.5 0.2;0 0 1];
colororderDash = [0.2 0.7 0.2;0.2 0.2 0.7;0.2 0.5 0.2;0 0 1];
set(h,'colororder',colororder) % verify so that legends becomes correct!!!

if 1 % Parallel electric field 1
    hca = h(isub); isub = isub + 1; 
    irf_plot(hca,{Epar3,Epar4},'comp')
    ylabel(hca,'E_{||} [mV/m]')
    set(hca,'colororder',colororder)
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end
%
if 1 % Perpendicular electric field 1
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,{Eper3,Eper4},'comp')
    ylabel(hca,'E_{\perp} [mV/m]')
    set(hca,'colororder',colororder)
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end
if 0 % Parallel magnetic field 1
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,BparAC3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,BparAC4,'color','b'); hold off;
    ylabel(hca,'\delta B_{||} [nT]')
    
    set(hca,'colororder',colororder)
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end
if 0 % Parallel magnetic field, model!
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,BB3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,BB4,'color','b'); hold off;
    ylabel(hca,'B_{||,model} [nT]')
    
    set(hca,'colororder',colororder)
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end
if 1 % Parallel magnetic field 1, including model!
    hca = h(isub); isub = isub + 1;  
    %irf_plot(hca,{BparAC3,BparAC4},'comp');
    set(hca,'colororder',colororderLight)
    reduc='*0.75';
    irf_plot(hca,{irf_filt(irf_tappl(BB3,reduc),fmin,fmax,fs,order),irf_filt(irf_tappl(BB4,reduc),fmin,fmax,fs,order)},'comp','linestyle','-'); hold(hca,'on');    
    
    set(hca,'colororder',colororderDash)
    irf_plot(hca,{irf_filt(irf_tappl(BB3,reduc),fmin,fmax,fs,order),irf_filt(irf_tappl(BB4,reduc),fmin,fmax,fs,order)},'comp','linestyle','--'); hold(hca,'on');
    set(hca,'colororder',colororder)
    irf_plot(hca,{irf_filt(staBpar3,fmin,fmax,fs,order),irf_filt(staBpar4,fmin,fmax,fs,order)},'comp');
    ylabel(hca,'\delta B_{||} [nT]')   
    
    
    alllines=findall(hca,'Type','line');
    realThick = 0;
    if realThick
        linew=1;
        set(alllines(1),'linewidth',linew)
        set(alllines(2),'linewidth',linew)   
    else
        linew=1;
        set(alllines(3),'linewidth',linew)
        set(alllines(4),'linewidth',linew)   
        set(alllines(5),'linewidth',linew)
        set(alllines(6),'linewidth',linew)
    end
    irf_legend(hca,{'C3','C4'},[0.98 0.9])    
    ylabel(hca,'\delta B_{||}[nT] ')
    userdata = get(gcf,'userdata');
    ylim = get(gca,'ylim');
    irf_legend(hca,'dashed - model',[0.06 0.92])
    %text(hca,0.1,ylim(1)+0.8*diff(ylim),'-- model')
end


% Finishing touches to figure
if printTitle; title(h(1),['fig6_staffonly_' num2str(fmin) 'Hz']); end
%title(h(1),['dt = ' num2str(abs(dt*1000),'%.0f') ' ms , v = ' num2str(v,'%.0f') ' km/s']);
irf_zoom(h,'x',tzoom{1});
%irf_zoom(h(2),'x',tzoom{1});
irf_plot_axis_align;
irf_timeaxis(hca,'usefig');

% Labeling
abcde={'a)','b)','c)', 'd)', 'e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.02,0.92],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
colorLegend={[0.2 0.5 0.2],[0 0 1]};
for k=1:nPanels%(isub-1)
    irf_legend(h(k),[abcde(k)],legloc{1},'color',colors{1})
    grid(h(k),'off')
    set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    %irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end

% Fix ylims for some axes
irf_zoom(h,'y')
%cn.print(['fig6_staffonly_' num2str(fmin) 'Hz'])
%% Make spectra of the noise and model holes for the shorter time interval
%staBpar3
%staBpar4
wds = 2.^(1:14);
wdstop = size(irf_tlim(staBpar3,tzoom{1}),1);
window = wds(find(wds<wdstop,1,'last'));
B3model = irf_resamp(BB3,staBpar3);
B4model = irf_resamp(BB4,staBpar4);
SR3a = irf_powerfft(irf_tlim(staBpar3,tzoom{1}),window,cn.f(staBpar3),80);
SR4a = irf_powerfft(irf_tlim(staBpar4,tzoom{1}),window,cn.f(staBpar3),80);
SR3b = irf_powerfft(irf_tlim(B3model,tzoom{1}),window,cn.f(B3model),80);
SR4b = irf_powerfft(irf_tlim(B4model,tzoom{1}),window,cn.f(B4model),80);
SR3.p_label = 'nT^2/Hz'; SR4.p_label = 'nT^2/Hz';
SR3a.pav = squeeze(nanmean(SR3a.p{1},1));
SR4a.pav = squeeze(nanmean(SR4a.p{1},1));
SR3b.pav = squeeze(nanmean(SR3b.p{1},1));
SR4b.pav = squeeze(nanmean(SR4b.p{1},1));
%%
h = axes; hold(h,'on')
plot(h,SR3a.f, SR3a.pav,'r');    
plot(h,SR4a.f, SR4a.pav,'g');    
plot(h,SR3b.f, SR3b.pav,'b');    
plot(h,SR4b.f, SR4b.pav,'k');  
set(h,'yscale','log','xscale','log')
legends = {'B3_{||,staff}','B4_{||,staff}','B3_{||,model}','B4_{||,model}'};
legend(h,legends)