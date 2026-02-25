% LOAD_LOCAL_DATA
dirData = '/Users/Cecilia/Data/MMS/2015Oct16/';

sc = 1:4;

% Spacecraft potential
c_eval('dobj_P?a = dataobj([dirData ''mms?/'' ''mms?_edp_brst_l2_scpot_20151016102944_v0.2.0.cdf'']);',sc)
c_eval('P?a = mms.variable2ts(get_variable(dobj_P?a,''mms?_edp_scpot''));',sc);
c_eval('dobj_P?b = dataobj([dirData ''mms?/'' ''mms?_edp_brst_l2_scpot_20151016103654_v0.2.0.cdf'']);',sc)
c_eval('P?b = mms.variable2ts(get_variable(dobj_P?b,''mms?_edp_scpot''));',sc);
%irf_plot({P1a,P1b},'comp')

% Burst electric fields
c_eval('dobj_E? = dataobj([dirData ''mms?/'' ''mms?_edp_brst_ql_dce2d_20151016102944_v0.4.0.cdf'']);',sc)
c_eval('dslE? = mms.variable2ts(get_variable(dobj_E?,''mms?_edp_dce_xyz_dsl''));',sc);
c_eval('dobj_P?b = dataobj([dirData ''mms?/'' ''mms?_edp_brst_l2_scpot_20151016103654_v0.2.0.cdf'']);',sc)
c_eval('P?b = mms.variable2ts(get_variable(dobj_P?b,''mms?_edp_scpot''));',sc);

%%

ts = mms.variable2ts(get_variable(dobj,dataVar{ii}));
ts = mms.variable2ts(get_variable(dobj,dataVar{ii}));
ts = mms.variable2ts(get_variable(dobj,dataVar{ii}));
ts = mms.variable2ts(get_variable(dobj,dataVar{ii}));

