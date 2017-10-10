function [phiV3 phiV4]= cn_phi_plot2(E3,E4,E3AC,E4AC,time,vhat,v,Te,t0)
% v in GSE
% E in ISR2
str='Yuris way';
% Getting the right velocity that correpsonds to X-ISR2 and Y-ISR2
[vx vy]=cn_v_isr2(time,v,vhat,2);
v3vec=ones(size(E3,1),2);
v3vec(:,1)=vx;
v3vec(:,2)=vy;
v4vec=ones(size(E4,1),2);
v4vec(:,1)=vx;
v4vec(:,2)=vy;

v=abs(v);
c_eval('E?_toint=irf_abs(E?AC);',3:4);
c_eval('phi?t=irf_integrate(E?_toint(:,1:3));',3:4);
c_eval('phiV?=[phi?t(:,1) -phi?t(:,2:3).*v?vec*6.2*1.6*0.1];',3:4);
c_eval('phi?xy=[phi?t(:,1) -phi?t(:,2:3).*v?vec*6.2*1.6*0.1/Te];',3:4);
c_eval('phi?=[phi?xy(:,1) sum(phi?xy(:,[2 3]),2)];',3:4);

c_eval('E?shift=E?;',4);
c_eval('E?shift(:,1)=E?shift(:,1)+t0;',4);

figure('name',['Potential (Yuri)'])
h=irf_plot(3);
isub=1;

lb=[0.5 0.8 1]; % light blue
lb=[0.7 0.8 1]; % light blue
g=[0.1 0.7 0]; % green

if 1 % Electric field
    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3(:,[1 2]),'color',g); hold(hca,'on');
    irf_plot(hca,E4(:,[1 2]),'b'); hold(hca,'on');
    irf_plot(hca,E4shift(:,[1 2]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 2]),'b'); hold(hca,'on');     
    irf_plot(hca,E3(:,[1 2]),'color',g); hold(hca,'on');
    ylabel(hca,'E_{x}[mV/m]');
    set(hca,'ColorOrder',[g;[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       
    irf_zoom(hca,'y');
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3(:,[1 3]),'color',g); hold(hca,'on');
    irf_plot(hca,E4(:,[1 3]),'b'); hold(hca,'on');
    irf_plot(hca,E4shift(:,[1 3]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 3]),'b'); hold(hca,'on'); 
    irf_plot(hca,E3(:,[1 3]),'color',g); hold(hca,'on'); 
    ylabel(hca,'E_{y}[mV/m]');
    set(hca,'ColorOrder',[g;[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
    irf_zoom(hca,'y');
end
if 1 % Potential
    hca=h(isub);isub=isub+1;
    irf_plot(hca,phi3(:,[1 2]),'color',g); hold(hca,'on');
    irf_plot(hca,phi4(:,[1 2]),'b'); hold(hca,'on');
    ylabel(hca,['e\phi/T_e  (',str,')']);
    set(hca,'ColorOrder',[g;[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]);
    set(hca,'ColorOrder',[0 0 0]);  
    irf_zoom(hca,'y');
end

tit=title(h(1),['Electric field and normalized electrostatic potential e\phi/T_ (ISR2)  v_{ph}=',num2str(v,'%.0f'),' km/s']);
%titpos=get(tit,'position');
%titpos(2)=titpos(2)+10.07;
%set(tit,'position',titpos)
irf_zoom(h,'x',[E3(1,1) E3(end,1)]);
%set(h,'XLimMode','manual')

%pos=get(h(1),'position')
%pos(4)=0.0001;
%pos(2)=pos(2)-0.008;
%set(h(1),'position',pos)
if 0
    
    xl=get(hca,'XLim');
    tcenter=mean(xl);
    distance=v*diff(xl)/2;
    logd=log10(distance);
    if logd>round(logd), 
        dx=10^(round(logd))/2;
    else
        dx=10^(round(logd))/5;
    end
    xticks=[-30:30]*dx/v/5+tcenter;
    xticklabels=cell(size(-30:30));
    xticklabels0=cell(size(-30:30));
    for j=-30:30, xticklabels{j+31}=' ';end
    for j=-6:6, xticklabels{j*5+31}=num2str(j*dx);end
    xticklabels;
    xticks;
    xticks0=get(h(1),'xtick');
    xmax=max(xticks0)+xticks0(end)-xticks0(end-1);
    index1=find(xticks>0,1,'first');
    index2=find(xticks<xmax,1,'last');

    xticklabels=xticklabels(index1:index2);
    xticks=xticks(index1:index2);
    midtick=xticks(ceil(end/2))
    midtick0=xticks0(ceil(end/2))
    xticklabels(ceil(end/2))
    max(xticks)
    xticks=xticks/midtick/2+(midtick0-midtick)
    2*(midtick0-midtick);
    if 1

    pos1=get(h(1),'position');
    pos1(3)=pos1(3);
    pos1(1)=pos1(1);
    pos1(2)=pos1(2);
    pos2=get(h(2),'position');

    ax2 = axes('Position',pos1,...%get(h(1),'Position'),...
    'XAxisLocation','top','YAxisLocation','right','Color','none',...
    'XColor','k','YColor','k','xtick',xticks*1/max(xticks),'xticklabel',xticklabels,...
    'ytick',get(h(1),'ytick'),'yticklabel',[]);
    xlabel(ax2,'km')
    end
end
if 1
    xticks0=get(h(1),'xtick');
    midindex=ceil(length(xticks0)/2);
    xticklabels=cell(length(xticks0));
    for l=1:length(xticks0); xticklabels{l}=' '; end 
    if mod(length(xticks0-1),4)==0; si=1; else si=2; end
    
    for l=si:2:length(xticks0)
        xticklabels{l}= num2str((xticks0(l)-xticks0(midindex))*v,'%.0f');
    end
    set(h(1),'xaxislocation','top','color','white','xtick',xticks0,'xticklabel',xticklabels);    
    xlabel(h(1),'km');
end
%set(h(1),'xaxislocation','top','color','none','xtick',xticks,'xticklabel',xticklabels);
%xlabel(h(1),'km');
%xgrid=get(hca,'xgrid');
%set(h(1),'xgrid',xgrid);
%set(h(2),'XLimMode','manual')
%set(h,'YLimMode','auto');
set(gcf,'PaperPositionMode','auto'); 
%irf_plot_axis_align
%irf_timeaxis(h,'usefig')
% eval(['print -dpdf ',t1str,'_',t2str,'_E_phi.pdf']);

