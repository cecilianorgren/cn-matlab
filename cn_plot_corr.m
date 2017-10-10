function cn_plot_corr(E3,E4,t0,str)
c_eval('E?shift=E?;',3:4);
c_eval('E?shift(:,1)=E?shift(:,1)-t0;',3:4);

figure('name',['E-shift: ',str])
h=irf_plot(4);
isub=1;

if 1 % Electric field
    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3(:,[1 2]),'g'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 2]),'b'); hold(hca,'on');
    ylabel(hca,'E_{x}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3shift(:,[1 2]),'g'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 2]),'b'); hold(hca,'on');
    ylabel(hca,'E_{x} (shifted) [mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end
if 1 % Electric field
    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3(:,[1 3]),'g'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 3]),'b'); hold(hca,'on');
    ylabel(hca,'E_{y}[mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,E3shift(:,[1 3]),'g'); hold(hca,'on');
    irf_plot(hca,E4(:,[1 3]),'b'); hold(hca,'on');
    ylabel(hca,'E_{y} (shifted) [mV/m]');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
    irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end


title(h(1),['Shifted electric field         t_{shift}=',num2str(t0),'   (',str,')'])
irf_zoom(h,'x',[E3(1,1) E3(end,1)]);