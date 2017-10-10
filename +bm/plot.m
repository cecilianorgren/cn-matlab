
tintDL = tailBMC3.TimeInterval(110,:);
tintStr = irf_time(tintDL(1),'yyyymmdd');
loadPath = ['/Users/Cecilia/Data/BM/' tintStr];

% load data
sc = 3;
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sc);
c_eval('diB?fgm=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sc);
%%
minvarE3 = irf_lmn(diE3,ud.v1,ud.v2,ud.v3);
minvarB3 = irf_lmn(diB3fgm,ud.v1,ud.v2,ud.v3);
%%
h = irf_plot(4);
isub=1;

hca=h(isub); isub=isub+1;
irf_plot(hca,diB3fgm)
ylabel(hca,'B_{ISR2}')
irf_legend(hca,{'x','y','z'},[0.95 0.9])

hca=h(isub); isub=isub+1;
irf_plot(hca,diE3)
ylabel(hca,'E_{ISR2}')
irf_legend(hca,{'x','y','z'},[0.95 0.9])

hca=h(isub); isub=isub+1;
irf_plot(hca,minvarB3)
ylabel(hca,'B_{MV}')
irf_legend(hca,{'x','y','z'},[0.95 0.9])

hca=h(isub); isub=isub+1;
irf_plot(hca,minvarE3)
ylabel(hca,'E_{MV}')
irf_legend(hca,{'x','y','z'},[0.95 0.9])

irf_zoom(h,'x',tintDL)
