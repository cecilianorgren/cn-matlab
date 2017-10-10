%% Reducing time series
tone=[2007 08 31 10 18 30 00];
ttwo=[2007 08 31 10 19 18 00];
B3zoom=cn_toepoch(tone,ttwo,absB3);
B4zoom=cn_toepoch(tone,ttwo,absB4);
E3zoom=cn_toepoch(tone,ttwo,diE3);
E4zoom=cn_toepoch(tone,ttwo,diE4);

%% Separating DC and AC E-field
dt=E3zoom(2,1)-E3zoom(1,1);
fs=1/dt;
flow=4;
fhigh=180;
c_eval('E?AC=irf_filt(E?zoom(:,1:3),flow,180,fs,3);',3:4);
c_eval('E?DC=irf_filt(E?zoom(:,1:3),0.000000001,flow,fs,3);',3:4);
c_eval('E?DC2=[E?zoom(:,1) E?zoom(:,[2 3])-E?AC(:,[2 3])];',3:4);



figure
h=irf_plot(6);
isub=1;
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);       

        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3AC(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4AC(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,AC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3AC(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4AC(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,AC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]); 
        
        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC2(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC2(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,DC}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);       

        hca=h(isub);isub=isub+1;
        irf_plot(hca,E3DC2(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4DC2(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X,DC2}[mV/m]');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);  
        

        title(h(1),['DC/AC E-field (ISR2)      cut at ',num2str(flow),'Hz'])
        irf_zoom(h,'x',[cn_toepoch(tone) cn_toepoch(ttwo)])
    irf_plot_axis_align
    add_timeaxis(hca,'usefig');
    
    
%%
E0=E3DC;
E0(:,[2 3])=0;
c_4_v_gui(E0,E0,E3DC,E4DC,2)
%%
E0=E3DC;
E0(:,[2 3])=0;
c_4_v_gui(E0,E0,E3AC,E4AC,2)
%%
    if 0 % Ex Ey filtered (2 panels)               
        c_eval('E?_filt=irf_filt(diE?zoom(:,1:3),flow,180,fs,3);',3:4);
        hca=h(isub); isub=isub+1;
        irf_plot(hca,E3_filt(:,[1 2]),'g'); hold(hca,'on');
        irf_plot(hca,E4_filt(:,[1 2]),'b'); hold(hca,'on');
        ylabel(hca,'E_{X}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]); 
        set(hca,'ColorOrder',[0 0 0]);       
        irf_legend(hca,{['f_{cut}=',num2str(flow),'Hz']},[0.02 0.90]); 
   
        hca=h(isub); isub=isub+1;
        irf_plot(hca,E3_filt(:,[1 3]),'g'); hold(hca,'on');
        irf_plot(hca,E4_filt(:,[1 3]),'b'); hold(hca,'on');
        ylabel(hca,'E_{Y}[mV/m]\newline ISR2');
        set(hca,'ColorOrder',[[0 1 0];[0 0 1]]);
        irf_legend(hca,{'C3','C4'},[0.02 0.05]);  
        set(hca,'ColorOrder',[0 0 0]);       
        irf_legend(hca,{['f_{cut}=',num2str(flow),'Hz']},[0.02 0.90]); 
    end