% examine electric field products
tintOverview = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tint = tintOverview;
%% Loading data
% ic = ;
ic = 1:4;
c_eval('dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint);',1:4);
c_eval('dslE?brst=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('dslE?3Dbrst=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('dslE?brst_l1b=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('dslE?3Dbrst_l1b=mms.db_get_ts(''mms?_edp_brst_l1b_dce'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('[EdB?,elang?] = irf_edb(irf.ts_vec_xy(dslE?brst.time,dslE?brst.data(:,1:2)),dmpaB?brst); EdB?.units = ''mV/m'';')
c_eval('[EparB?,elang?] = irf_edb(irf.ts_vec_xy(dslE?brst.time,dslE?brst.data(:,1:2)),dmpaB?brst.resample(dslE?brst.time),20,''Epar''); EparB?.units = ''mV/m'';')
c_eval('EparB?.data(abs(elang?.data)>20,3)=NaN;')
c_eval('EdB?.data(abs(elang?.data)<20,3)=NaN;')
%%
ic = 4;
h = irf_plot(6);

hca = irf_panel('B');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('brst Ex');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?3Dbrst.tlim(tint).x,EdB?.x},''comp'');',ic)
hca.YLabel.String = {'E_x','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'dce2d','dce','E dot B = 0'},[0.95 0.95]);

hca = irf_panel('brst Ey');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).y,dslE?3Dbrst.tlim(tint).y,EdB?.y},''comp'');',ic)
hca.YLabel.String = {'E_y','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'dce2d','dce','E dot B = 0'},[0.95 0.95]);

hca = irf_panel('brst Ez');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).z,dslE?3Dbrst.tlim(tint).z,EdB?.z,EparB?.z},''comp'');',ic)
hca.YLabel.String = {'E_z','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'dce2d','dce','E dot B = 0','E dot B ~= 0'},[0.95 0.95]);

if 0
  hca = irf_panel('brst Eabs');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?brst.tlim(tint).abs,dslE?3Dbrst.tlim(tint).abs,EdB?.abs,EparB?.abs},''comp'');',ic)
  hca.YLabel.String = {'|E|','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'dce2d','dce','E dot B = 0','E dot B ~= 0'},[0.95 0.95]);
end

hca = irf_panel('brst Ez diff');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).z-dslE?3Dbrst.tlim(tint).z,dslE?brst.tlim(tint).z-EdB?.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E_{z}','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_{z,dce2d}-E_{z,dce}','E_{z,dce2d}-E_{z,E dot B = 0}'},[0.95 0.95]);

hca = irf_panel('elevation angle');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{elang1,elang2,elang3,elang4},'comp')
hca.YLabel.String = {'El. angle','(deg)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);
hca.YTick = [-90 -20 0 20 90];

h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_axis_align
irf_zoom(h,'y')
irf_zoom(h,'x',tint)