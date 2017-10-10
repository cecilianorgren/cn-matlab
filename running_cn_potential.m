gseEac=irf_filt(cn_toepoch(t1,t2,diE3),flow,0,fs,3);
potential=cn_potential(E3AC,'dsi',vd*vchoice,'gse');
potential2=cn_potential(gseEac,'gse',vd*vchoice,'gse');
figure;irf_plot({potential,potential2,phi3,phiV3},'comp')
legend('dsi','gse','yuri','andris')
phiTeav=irf_multiply(1/Teav,potential,0.5,potential,0.5);
phiTeav2=[potential(:,1) potential(:,2:end)/Teav];

figure;irf_plot({phiTeav,phiTeav2},'comp')

%%
c_eval('e?=cn_toepoch(t1,t2,diE?);',3:4)
c_eval('b?=cn_toepoch(t1,t2,gseB?);',3:4)
c_eval('epar?=irf_edb(e?,b?,90,''Epar'');',3:4)
c_eval('epar?=c_coord_trans(''dsi'',''gse'',epar?,''cl_id'',?);',3:4)
c_eval('eac?=irf_filt(epar?,5,0,450,3);',3:4)
%%
tguess=0.009;
cn_plot_corr(e3,e4,tguess,'gissat');
vguess=-zdist/tguess;
%%
c_eval('b?abs=irf_abs(b?);',3:4)
c_eval('khat?=[b?abs(:,1) b?abs(:,2:4)./[b?abs(:,5) b?abs(:,5) b?abs(:,5) ]];',3:4)
c_eval('pot?=cn_potential(eac?,''gse'',khat?*vguess,''gse'');',3:4)
%figure; irf_plot({epar3,epar4},'comp')
%figure; irf_plot({pot3,pot4},'comp')

%%
figure('name','epar');
h=irf_plot(5);
isub=1;
c_eval('epar?=eac?;',3:4)
if 1 % Electric field
    lb=[0.5 0.8 1]; % light blue
    lb=[0.6 0.7 1]; % light blue
    gr=[0 0.8 0]; % green
    c_eval('epar?=irf_abs(epar?);',3:4);
    c_eval('epars?=epar?; epars?(:,1)=epars?(:,1)+tguess;',4);
    
    hca=h(isub);isub=isub+1;
    irf_plot(hca,epar3(:,[1 2]),'color',gr); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 2]),'b'); hold(hca,'on');
    irf_plot(hca,epars4(:,[1 2]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 2]),'b'); hold(hca,'on');  
    irf_plot(hca,epar3(:,[1 2]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{x}[mV/m] GSE');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);       

    hca=h(isub);isub=isub+1;
    irf_plot(hca,epar3(:,[1 3]),'color',gr); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 3]),'b'); hold(hca,'on');
    irf_plot(hca,epars4(:,[1 3]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 3]),'b'); hold(hca,'on');  
    irf_plot(hca,epar3(:,[1 3]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{y}[mV/m] GSE');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);     

    hca=h(isub);isub=isub+1;
    irf_plot(hca,epar3(:,[1 4]),'color',gr); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 4]),'b'); hold(hca,'on');
    irf_plot(hca,epars4(:,[1 4]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 4]),'b'); hold(hca,'on');  
    irf_plot(hca,epar3(:,[1 4]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{z}[mV/m] GSE');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end
if 1
    c_eval('eabs?=cn_m_trans(epar?,M,1);',3:4)
    c_eval('eabss?=epar?; eabss?(:,1)=eabss?(:,1)+tguess;',4);
    hca=h(isub);isub=isub+1;
    irf_plot(hca,eabs3(:,[1 4]),'color',gr); hold(hca,'on');
    irf_plot(hca,eabs4(:,[1 4]),'b'); hold(hca,'on');
    irf_plot(hca,eabss4(:,[1 4]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,eabs4(:,[1 4]),'b'); hold(hca,'on');  
    irf_plot(hca,eabs3(:,[1 4]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{||}[mV/m] GSE');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end
if 0
    hca=h(isub);isub=isub+1;
    irf_plot(hca,epar3(:,[1 5]),'color',gr); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 5]),'b'); hold(hca,'on');
    irf_plot(hca,epars4(:,[1 5]),'color',lb,'linestyle','--'); hold(hca,'on');
    irf_plot(hca,epar4(:,[1 5]),'b'); hold(hca,'on');  
    irf_plot(hca,epar3(:,[1 5]),'color',gr); hold(hca,'on');
    irf_zoom(hca,'y');
    
    ylabel(hca,'E_{||}[mV/m] GSE');
    set(hca,'ColorOrder',[[0 1 0];[0 0 1];lb]);
    irf_legend(hca,{'C3','C4','C4_{shifted}'},[0.02 0.05]); 
    set(hca,'ColorOrder',[0 0 0]);  
end

if 1 % Potential and dB
    hca=h(isub); isub=isub+1;
    c_eval('pote?=[pot?(:,1) pot?(:,2)/Teav];',3:4)
 
    
    irf_plot(hca,pote3,'color',[0.8 0.2 0.2]); hold on;    
    irf_plot(hca,pote4,'color',[0.8 0.8 0.2]); hold on;
    irf_zoom('x',[cn_toepoch(t1) cn_toepoch(t2)]);
    set(hca,'ColorOrder',[[0.8 0.2 0.2];[0.8 0.8 0.2]]);
    irf_legend({'C3','C4'},[0.02 0.05])
    ylabel('e\phi/T_e')
end

title(h(1),['Parallel electric field and potential    (T_e=',num2str(Teav,'%.0f'),' eV   v=',num2str(vguess,'%.0f'),' km/s)'])
irf_zoom(h,'x',[epar3(1,1) epar3(end,1)])
