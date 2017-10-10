%% Load mms data from local files
cd /Users/Cecilia/Data/MMS/20150815/
c_eval('mms?_edp_brst_ql_dce=dataobj(''mms?_edp_brst_ql_dce*.cdf'');',3)

c_eval(['E? = mms.variable2ts(get_variable(mms?_edp_brst_ql_dce,''mms?_edp_dce_xyz_dsl''));'...
        'E?.userData.LABLAXIS=''E''; E?.coordinateSystem=''dsl'';'],3)

%tint = irf.tint('2015-08-15T13:10:00Z/2015-08-15T13:50:00Z');
%c_eval('E? = mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint);',3);
%%
sc = 3;
c_eval('mms?_edp_brst_l1b_dce=dataobj(''mms?_edp_brst_l1b_dce*.cdf'');',sc)
e12 = mms.variable2ts(get_variable(mms3_edp_brst_l1b_dce,'mms3_edp_e12_enable'));
c_eval(['e12 = mms.variable2ts(get_variable(mms?_edp_brst_l1b_dce,''mms?_edp_e12_enable''));'...
  'e12.userData.LABLAXIS=''E_{12}''; e12.coordinateSystem=''dsl'';'],sc)

%%
c_eval('mms?_dfg_srvy_ql=dataobj(''mms?_dfg_srvy_ql*.cdf'');')
c_eval(['B? =mms.variable2ts(get_variable(mms?_dfg_srvy_ql,''mms?_dfg_srvy_dmpa''));'...
  'B?.userData.LABLAXIS=''B''; B?.coordinateSystem=''dmpa'';'])
%% sc potential
sc=1:4;
c_eval('mms?_edp_brst_l2_scpot=dataobj(''mms?_edp_brst_l2_scpot_20150815130004_v0.1.0.cdf'');',sc)
c_eval('scp? = mms.variable2ts(get_variable(mms?_edp_brst_l2_scpot,''mms?_edp_scpot''));',sc);

%% magnetic field?
%sc=3;
%c_eval('mms?_dfg_slow_l1a=dataobj(''mms?_dfg_slow_l1a_20150815_v0.7.3.cdf'');',sc)
%c_eval('dfg? = mms.variable2ts(get_variable(mms?_dfg_slow_l1a,''mms?_dfg_123''));',sc);
sc=1:4;
c_eval('mms?_dfg_brst_ql=dataobj(''mms?_dfg_brst_ql_20150815125500_v0.0.0.cdf'');',sc)
c_eval('R? = mms.variable2ts(get_variable(mms?_dfg_brst_ql,''mms?_ql_pos_gsm''));',sc);

c_eval('R?=mms?_dfg_brst_ql.data.mms?_ql_pos_gsm.data;',sc)
%%
h=irf_plot({E3,scp3});
h(1).YLim = 40*[-1 1];
%%
h = irf_plot(2,'newfigure');
hca = irf_panel('E');
irf_plot(hca,E3)
hca.YLim = 40*[-1 1];
hca = irf_panel('Vsc');
irf_plot(hca,{scp1,scp2,scp3,scp4},'comp')


%% Load data
tint = irf.tint('2015-08-15T13:00:00Z/2015-08-15T13:04:00Z');

load /data/mms/irfu/mmsR.mat
epoTmp = EpochTT(R.time);
gsmR1 = [epoTmp.epochUnix R.gsmR1];

scId = 3;
c_eval([...
    %'E? = mms.db_get_ts (''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint);'...
    'P? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_psp'',tint);'...
    %'B? = mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_gsm_dmpa'',tint);'...
    ],scId)

%% Plot
% define Cluster colors
mmsColors=[0 0 0; 1 0 0 ; 0 0.5 0 ; 0 0 1];
% Official MMS colors
%mmsColors=[0 0 0; .8 .4 0 ; 0 0.6 0.5 ; 0.35 0.7 .9];

h = irf_plot(3,'newfigure');

hca = irf_panel('B'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,B3)
ylabel(hca,'B3 [nT]')

hca = irf_panel('E'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,E3)
ylabel(hca,'E3 [mV/m]')

hca = irf_panel('P'); set(hca,'ColorOrder',mmsColors)
irf_plot(hca,P3)
ylabel(hca,'P3 [V]')
%%

irf_zoom(h,'x',tint)
irf_plot_axis_align(h)
add_position(h(end),gsmR1)
xlabel(h(end),'')
title(h(1),tint.start.utc)

%% Load burst data from local files
scId = 3;
fileTimeB = '20150815125500';
c_eval('dobjB?=dataobj([''mms?_dfg_brst_ql_20150815125500_v0.0.0.cdf'']);',scId)
c_eval(['B?=mms.variable2ts(get_variable(dobjB?,''mms?_dfg_brst_gsm_dmpa''));'...
        'B?.userData.LABLAXIS=''B''; B?.coordinateSystem=''dmpa'';'],scId)
    
c_eval('dobjP?=dataobj(''mms?_edp_brst_l2_scpot_20150815132504_v0.1.0.cdf'');',scId)
c_eval(['P?=mms.variable2ts(get_variable(dobjP?,''mms?_edp_dcv''));'...
        'P?.userData.LABLAXIS=''P'';'],scId)
        
c_eval('dobjE?=dataobj(''mms?_edp_brst_ql_dce2d_20150815132504_v0.1.0.cdf'');',scId)
c_eval(['E?=mms.variable2ts(get_variable(dobjE?,''mms?_edp_dce_xyz_dsl''));'...
        'E?.userData.LABLAXIS=''E''; E?.coordinateSystem=''DSL'';'],scId)
    
c_eval('dobjB?scm=dataobj(''mms?_scm_brst_l1a_scb_20150815132504_v0.8.0.cdf'');',scId)
c_eval(['B?scm=mms.variable2ts(get_variable(dobjB?scm,''mms?_scm_scb_123''));'...
        'B?scm.userData.LABLAXIS=''B (SCM)''; B?scm.coordinateSystem=''?'';'],scId)
    
c_eval('dobjEbot?=dataobj(''mms?_feeps_srvy_l1a_electron-bottom_20150815000000_v2.2.3.cdf'');',scId)
c_eval(['ecounts?=mms.variable2ts(get_variable(dobjEbot?,''mms?_epd_feeps_survey_live_time_counts''));'...
        'ecounts?.userData.LABLAXIS=''electron counts''; ecounts?.coordinateSystem=''sensorID 1'';'],scId)    