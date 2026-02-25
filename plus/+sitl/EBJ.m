tint = irf.tint('2015-12-01T00:00:00.00Z/2015-12-01T14:00:00.00Z');

% load data 
ic = 1:4;

c_eval('dslE?=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);

c_eval('dmpaB?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);

c_eval('gseR?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_ql_pos_gse'',tint);',ic);

c_eval('ne?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);

[j,divB,B,jxB,divTshear,divPb] = c_4_j('gseR?','dmpaB?');
j = irf.ts_vec_xyz(j.time,j.data*1e9);
j.units = 'nAm^{-2}';

%%
h = irf_plot(7); 
isub = 1;

hca = irf_panel('Bx');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.y.tlim(tint),dmpaB2.y.tlim(tint),dmpaB3.y.tlim(tint),dmpaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);

hca = irf_panel('Bz');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z}','(nT)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'sc1','sc2','sc3','sc4'},[1.01 0.9]);


hca = irf_panel('J');
set(hca,'ColorOrder',mms_colors('xyza'))
irf_plot(hca,{j.x,j.y,j.z},'comp');
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J','[nA m^{-2}]'},'Interpreter','tex');
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'J_x','J_y','J_z'},[0.95 0.95]);

hca = irf_panel('E');
icpl = 3;
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?.tlim(tint).x,dslE?.tlim(tint).y,dslE?.tlim(tint).z},''comp'');',icpl)
hca.YLabel.String = {irf_ssub('E mms?',icpl),'(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[0.95 0.95]);

hca = irf_panel('ne');
irf_plot(ne1,'comp')
irf_zoom(h,'x',tint)


