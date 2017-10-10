
t1=[2007 09 26 10 13 00 00];
t2=[2007 09 26 10 30 00 00];
t1=[2007 09 26 09 48 00 00];
t2=[2007 09 26 09 56 00 00];
%%
t1=[2007 09 2 14 30 00 00];
t2=[2007 09 2 14 40 00 00];
t=[2007 09 2 14 30 29 00];
tint=[cn_toepoch(t1) cn_toepoch(t2)];
%%

diE3z=cn_toepoch(t1,t2,diE3);
diE4z=cn_toepoch(t1,t2,diE4);

figure;h=irf_plot(5);
isub=1;
c_eval('diB?=c_coord_trans(''gsm'',''dsi'',gsmB?,''cl_id'',?);',3:4);
deltaB=irf_add(1,diB3,-1,diB4);


Ne_t=cn_toepoch(t,peaNe4);
    Ne=Ne_t(2)*1e6;
    mu0=4*pi*1e-7;
    e=1.6e-19; 
    B0_t=cn_toepoch(t,gsmB4);
    B0=B0_t(5);
    eVTe_t=cn_toepoch(t,eVTe3);
    eVTe=eVTe_t(2);
    scaling=(Ne*mu0*e/B0)*1e18;
    flh_freq=irf_plasma_calc(B0,Ne*1e-6,0,0,0,'Flh');
    deltaB4=irf_filt(gsmB4,0.5*flh_freq,0,450,3);
dphi=[deltaB4(:,1) deltaB4(:,[2])/scaling/eVTe];

B4DC=irf_add(1,gsmB4,-1,deltaB4);


if 0
    hca=h(isub);isub=isub+1;
irf_plot(h(1),diE3z(:,[1 3]),'k');hold on;
irf_plot(h(1),diE4z(:,[1 3]),'r');hold on;
irf_legend(h(1),{'C3','C4'},[0.02 0.04])
ylabel(h(1),'E_x [mV/m] ISR2')
end
if 1
    hca=h(isub);isub=isub+1;
irf_plot(hca,gsmB4);hold on;
irf_legend(hca,{'x','y','z','|B|'},[0.02 0.04])
ylabel(hca,'B [nT] GSM C4')
end
if 0
    hca=h(isub);isub=isub+1;
irf_plot(hca,B4DC);hold on;
irf_legend(hca,{'x','y','z','|B|'},[0.02 0.04])
ylabel(hca,'B_0 [nT] GSM C4')
end
if 1
    hca=h(isub);isub=isub+1;
irf_plot(hca,deltaB4(:,1:4));hold on;
irf_legend(hca,{'x','y','z',},[0.02 0.04])
ylabel(hca,'B_1 [nT] GSM C4')
end
if 1
    hca=h(isub);isub=isub+1;
irf_plot(hca,deltaB(:,1:4));hold on;
irf_legend(hca,{'x','y','z',},[0.02 0.04])
ylabel(hca,'\Delta B [nT] ISR2 C4')
end
if 1
    hca=h(isub);isub=isub+1;
    
irf_plot(hca,dphi);
irf_legend(hca,{'phi from B_{1,x}',},[0.02 0.04])
ylabel(hca,'e\phi/T_e')
end
if 0
        hca=h(isub);isub=isub+1;
irf_plot(hca,deltaB4(:,[1 2]),'k');
irf_legend(hca,{'C4'},[0.02 0.04])
ylabel(hca,'B_1 [nT]')
end
if 1
    hca=h(isub);isub=isub+1;
irf_plot(hca,diE3z(:,[1 3]),'k');hold(hca, 'on');
irf_plot(hca,diE4z(:,[1 3]),'r');hold(hca, 'on');
irf_legend(hca,{'C3','C4'},[0.02 0.04])
ylabel(hca,'E_y [mV/m] ISR2')
end


if 0
irf_plot(h(3),diE3z(:,[1 3]),'k');hold(hca, 'on');
irf_plot(h(3),diE4z(:,[1 3]),'r');hold(hca, 'on');
irf_legend(h(3),{'C3','C4'},[0.02 0.04])
ylabel(h(3),'E_y [mV/m] ISR2')
end
title(h(1),['n_e=',num2str(Ne_t(2),'%.2f'),' cm^{-3}  T_e=',num2str(eVTe,'%.0f'),' eV   \omega_{LH}=',num2str(flh_freq,'%.0f'),' Hz   \rho_e=1.5 km'])
irf_zoom(h,'x',tint)