tint = irf.tint('2015-10-08T06:00:00.000Z/2015-10-08T16:30:00.000Z');
sc=1:4;

%% Load some survey data
%c_eval('[scm?,dobj_scm?] = mms.cn_get_ts(''mms?_scm_fast_l2_scf'',''mms?_scm_scf_gse'',tint(1));',sc);
c_eval('[edp?,dobj_edp?] = mms.cn_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint(1));',sc);
%c_eval('[fpi?,dobj_fpi?] = cn_get_ts(''mms1_fpi_fast_sitl'',''mms?_edp_dce_xyz_dsl'',tint(1));',sc);
%c_eval('[pa_mid?,dobj_fpi?] = cn_get_ts(''mms1_fpi_fast_sitl'',''mms?_fpi_ePitchAngDist_midEn'',tint(1));',sc);
c_eval('[dfg?,dobj_dfg?] = cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint(1));',sc)

%% Plot and compare

%tint_zoom = irf.tint('2015-10-16T06:00:00',8*60*60);
h = irf_plot({dfg1,edp1,dfg2,edp2,dfg3,edp3,dfg4,edp4});
irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_pl_mark(h,tint(1).epochUnix')


%%
%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
h = irf_plot(6);
hca = irf_panel('Bx');
irf_plot(hca,{dfg1.x,dfg2.x,dfg3.x,dfg4.x},'comp');
hca.YLabel.String = 'B_x [nT]';
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('By');
irf_plot(hca,{dfg1.y,dfg2.y,dfg3.y,dfg4.y},'comp');
hca.YLabel.String = 'B_y [nT]';
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz');
irf_plot(hca,{dfg1.z,dfg2.z,dfg3.z,dfg4.z},'comp');
hca.YLabel.String = 'B_z [nT]';
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ex');
irf_plot(hca,{edp1.x,edp2.x,edp3.x,edp4.x},'comp');
hca.YLabel.String = 'E_x [nT]';
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ey');
irf_plot(hca,{edp1.y,edp2.y,edp3.y,edp4.y},'comp');
hca.YLabel.String = 'E_y [nT]';
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ez');
irf_plot(hca,{edp1.z,edp2.z,edp3.z,edp4.z},'comp');
hca.YLabel.String = 'E_z [nT]';
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);


irf_zoom(h,'x',tint_zoom)
irf_zoom(h,'y')
irf_pl_mark(h,tint(1).epochUnix')

% Roy Torberts requested intervals
tint1 = irf.tint('2015-10-16T13:00:00',8*60);
tint2 = irf.tint('2015-10-16T11:25:00',60*7);
irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')
