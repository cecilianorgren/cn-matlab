%% Adding multiple axes to irf_plot
h=irf_plot(2);
isub=1;
hca=h(isub); isub=isub+1;
irf_plot(hca,peaNe3); hold(hca,'on')
%ax1=hca;
%set(ax1,'XColor','r','YColor','r')
%ud=get(ax1
%%       
hca=h(isub);
ud = get(gcf,'userdata');
t_st_e = double(ud.t_start_epoch);
        
plot(hca,peaNe3(:,1)-t_st_e,peaNe3(:,2)); hold(hca,'on');
ax1=hca;
%set(ax1,'xticklabel',[],'xtick',[]);
ax2=axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k',...
           'xticklabel',[],'xtick',[]);
  
plot(ax2,(betaEH(:,1)-t_st_e), betaEH(:,2),'k','Parent',ax2)
%ylabel(hca,'n_e [cm^{-3}]');
%set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
set(ax1,'xticklabel',[],'xtick',[],...
    'yscale','log','yaxislocation','right')
irf_legend(hca,{'PEACE'},[0.02 0.2]);
%%       
h=irf_plot(2);isub=1;
hca=h(isub);
%ud = get(gcf,'userdata');
%t_st_e = double(ud.t_start_epoch);
        
irf_plot(hca,EP); 
set(hca,'box','off')
ax1=hca;
set(ax1,'ycolor','r')
%set(ax1,'xticklabel',[],'xtick',[]);
ax2=axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k',...
           'xticklabel',[],'xtick',[],...
           'yscale','log');
  
irf_plot(hca,IPNT,'k','Parent',ax2)
%ylabel(hca,'n_e [cm^{-3}]');
%set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);

%set(ax1,'xticklabel',[],'xtick',[],...
%    'yscale','log','yaxislocation','right')
irf_legend(hca,{'PEACE'},[0.02 0.2]);
%% Matlab example
figure;
x1 = [0:.1:40];
y1 = 4.*cos(x1)./(x1+2);
x2 = [1:.2:20];
y2 = x2.^2./x2.^3;

hl1 = line(x1,y1,'Color','r');
ax1 = gca;
set(ax1,'XColor','r','YColor','r')

ax2 = axes('Position',get(ax1,'Position'),...
           'XAxisLocation','top',...
           'YAxisLocation','right',...
           'Color','none',...
           'XColor','k','YColor','k');
hl2 = line(x2,y2,'Color','k','Parent',ax2);
%% irfnotes example
h=irf_plot(1);
% plot data
irf_plot(h(1),EP)
set(h(1),'box','off')
%
h(2) = axes('Position',get(h(1),'Position'));
irf_plot(h(2),IPNT,'r')
set(h(2),'XAxisLocation','top','xtick',[]); % remove 'xtick' if xticks required
set(h(2),'YAxisLocation','right','YScale','log');
set(h(2),'Color','none','box','off'); % color of axis
set(h(2),'XColor','r','YColor','r'); % color of axis lines and numbers

irf_timeaxis(h(2),'nolabels')

irf_legend(h(1),'data',[0.02 0.98],'color','k')
irf_legend(h(2),'sqrt(data)',[0.8 0.98],'color','r')

ylabel(h(1),'data')
ylabel(h(2),'sqrt(data)')