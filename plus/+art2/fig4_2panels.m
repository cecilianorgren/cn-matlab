% This script produces a figure with either one or two panels with
% wavetrains, that serve to illustrate the high probability that we observe
% the same holes on both sc.

% fig4: a wave train
cd /Users/Cecilia/Data/BM/20070831
[tint, quality, comments]=eh_tint;
highQuality = find(quality==2); % highest quality fields

%if ~exist('ind','var'); ind = 6; end
%
%tint = tint{highQuality(ind)};
tzoom{1} = toepoch([2007 08 31 10 17 45.3;2007 08 31 10 17 46.9])'; % ind=5, quality=2
tzoom{2} = toepoch([2007 08 31 10 17 51.2;2007 08 31 10 17 52.4])';
%tzoom{2} = tint{highQuality(11)};
%tzoom=tint;


% Make figure
set(0,'defaultLineLineWidth', 1);
fig = figure(24);
set(gcf,'defaultAxesFontSize',14);
set(gcf,'defaultTextFontSize',14);
set(gcf,'defaultAxesFontUnits','pixels');
set(gcf,'defaultTextFontUnits','pixels');
set(fig,'position',[10   684   591   195]);
nPanels = 2;
h(1) = subplot(2,1,1); set(h(1),'position',[0.1300    0.65    0.7750    0.33]);
h(2) = subplot(2,1,2); set(h(2),'position',[0.1300    0.2100    0.7750    0.33])

isub = 1;

if 1 % Parallel electric field
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,Epar3,'color',[0.2 0.5 0.2]); hold(hca,'on');
    irf_plot(hca,Epar4,'color','b'); hold off;
    ylabel(hca,'E_{||} [mV/m]')
    set(hca,'colororder',[1 1 0; 0 0 1])
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end

if 1 % Parallel electric field
    hca = h(isub); isub = isub + 1;  
    irf_plot(hca,Epar3,'color',[0.2 0.5 0.2]); hold on;
    irf_plot(hca,Epar4,'color','b'); hold off;
    ylabel(hca,'E_{||} [mV/m]')
    set(hca,'colororder',[1 1 0; 0 0 1])
    irf_legend(hca,{'C3','C4'},[0.98 0.9])
end

%xlabel(h(1),'')

% Finishing touches to figure
%title(h(1),['dt = ' num2str(abs(dt*1000),'%.0f') ' ms , v = ' num2str(v,'%.0f') ' km/s']);
irf_zoom(h(1),'x',tzoom{1});
irf_zoom(h(2),'x',tzoom{2});
irf_plot_axis_align;
irf_timeaxis(hca,'usefig');

% Labeling
abcde={'', 'b)', 'c)', 'd)', 'e)','f)','g)','h)','i)','j)'};
colors={[0 0 0],[0 0 0],[1 1 1],[0 0 0],[0 0 0]};
legloc={[0.02,0.92],[0.02,0.08],[0.02,0.08],[0.02,0.96],[0.02,0.96],[0.02,0.96],[0.02,0.96]};
colorLegend={[0.2 0.5 0.2],[0 0 1]};
for k=1:(isub-1)%nPanels
    irf_legend(h(k),[abcde(k)],legloc{1},'color',colors{1})
    grid(h(k),'off')
    set(h(k),'colororder',[colorLegend{1}; colorLegend{2}])
    irf_legend(h(k),{'C3','C4'},[0.98 0.9])
end

% Fix ylims for some axes
irf_zoom(h,'y')

