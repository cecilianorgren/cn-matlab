% plot of wavetrain with measured parallel magnetic field
% load data from art2.load_data

    
%% Make plot !!!
% This script produces a figure with either one or two panels with
% wavetrains, that serve to illustrate the high probability that we observe
% the same holes on both sc.
printTitle = 0;
% fig4: a wave train
cd /Users/Cecilia/Data/BM/20070831
tzoom{1} = toepoch([2007 08 31 10 17 45.3;2007 08 31 10 17 46.9])'; % ind=5, quality=2
tzoom{2} = toepoch([2007 08 31 10 17 51.2;2007 08 31 10 17 52.4])';
fmin = 1.2;
fmax = 0;
fs = cn.f(staBpar3);
order = 5;
if 1 % filter
    B3 = irf_filt(staBpar3,fmin,fmax,fs,order);
    B4 = irf_filt(staBpar4,fmin,fmax,fs,order);
else % detrend
    B3 = [staBpar3(:,1), detrend(staBpar3(:,2))];
    B4 = [staBpar4(:,1), detrend(staBpar4(:,2))];
end
    
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
if 1 % Parallel magnetic field 1
    hca = h(isub); isub = isub + 1;  
    set(hca,'colororder',colororder)
    irf_plot(hca,{B3,B4},'comp');     
    ylabel(hca,'\delta B_{||} [nT]')      
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
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
legloc={[0.02,0.95],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
colorLegend={[0.2 0.5 0.2],[0 0 1]};
for k=1:nPanels%(isub-1)
    irf_legend(h(k),[abcde(k)],legloc{1},'color',colors{1})
    grid(h(k),'off')
    %set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    %irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end

% Fix ylims for some axes
irf_zoom(h,'y')
%cn.print(['fig6_staffonly_' num2str(fmin) 'Hz'])
