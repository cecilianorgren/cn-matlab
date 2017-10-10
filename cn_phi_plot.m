function [phiV3 phiV4] = cn_phi_plot(E3,E4,E3AC,E4AC,time,vd,v,Te,tref,t1str,t2str,str,t0)
%% Getting the right velocity that correpsonds to X-ISR2 and Y-ISR2
[vx vy]=cn_v_isr2(time,v,vd,1)
v3vec=ones(size(E3,1),2);
v3vec(:,1)=vx;
v3vec(:,2)=vy;
v4vec=ones(size(E4,1),2);
v4vec(:,1)=vx;
v4vec(:,2)=vy;

v=abs(v);
c_eval('E?_toint=irf_abs(E?AC);',3:4);
c_eval('phi?t=irf_integrate(E?_toint(:,1:3));',3:4);
c_eval('phi?=[phi?t(:,1) -phi?t(:,2:3).*v?vec*6.2*1.6*0.1/Te];',3:4);
c_eval('phiV?=[phi?t(:,1) -phi?t(:,2:3).*v?vec];',3:4);
c_eval('E?shift=E?; E?shift(:,1)=E?shift(:,1)-t0;',3);
figure('name',['Potential (Andris)'])
h=irf_plot(4);
isub=1;

lb=[0.5 0.8 1]; % light blue
lb=[0.7 0.8 1]; % light blue
lb=[0.5 1 0.2]; % light blue
g=[0 0.8 0]; % green

if 1 % Electric field
    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3(:,[1 2]),'color',g); hold(hca,'on');
    irf_plot(hca,E4(:,[1 2]),'b'); hold(hca,'on');
    irf_plot(hca,E3shift(:,[1 2]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 2]),'b'); hold(hca,'on');  
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{x}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C3_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3(:,[1 3]),'color',g); hold(hca,'on');
    irf_plot(hca,E4(:,[1 3]),'b'); hold(hca,'on');
    irf_plot(hca,E3shift(:,[1 3]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 3]),'b'); hold(hca,'on');  
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{y}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C3_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end
if 1 % Potential
    hca=h(isub);isub=isub+1;
    irf_plot(hca,phi3(:,[1 2]),'color',g); hold(hca,'on');
    irf_plot(hca,phi4(:,[1 2]),'b'); hold(hca,'on');
    ylabel(hca,['e\phi_{x}/T_e  (',str,')']);
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]);  
    irf_zoom(hca,'y');

    hca=h(isub);isub=isub+1;
    irf_plot(hca,phi3(:,[1 3]),'color',g); hold(hca,'on');
    irf_plot(hca,phi4(:,[1 3]),'b'); hold(hca,'on');
    ylabel(hca,['e\phi_{y}/T_e  (',str,')']);
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]);  
    irf_zoom(hca,'y');
end

title(h(1),['Electric field and normalized electrostatic potential e\phi/T_ (ISR2)  v_{ph}=',num2str(v,'%.0f'),' km/s'])
irf_zoom(h,'x',[E3(1,1) E3(end,1)]);
set(h,'XLimMode','manual')
%irf_plot_axis_align
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
set(h(1),'xaxislocation','top','xtick',xticks,'xticklabel',xticklabels);
xlabel(h(1),'km');
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
set(h(2),'XLimMode','manual')
%set(h,'YLimMode','auto');
set(gcf,'PaperPositionMode','auto'); 
eval(['print -dpdf ',t1str,'_',t2str,'_E_phi.pdf']);
