tint = irf.tint('2015-10-16T10:00:00.00Z/2015-10-16T14:00:00.00Z');
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.60Z');
ic = 1;

%% Burst data
Tint = tint;
c_eval('datasetName = ''mms?_fpi_brst_l2_dis-dist'';',ic)
c_eval('pref = ''mms?_dis'';',ic)
Vr.tmmode = 'brst';
dist = mms.db_get_ts(datasetName,[pref '_dist_' Vr.tmmode],Tint);
energy0 = mms.db_get_variable(datasetName,[pref '_energy0_' Vr.tmmode],Tint);
energy1 = mms.db_get_variable(datasetName,[pref '_energy1_' Vr.tmmode],Tint);
phi = mms.db_get_ts(datasetName,[pref '_phi_' Vr.tmmode],Tint);
theta = mms.db_get_variable(datasetName,[pref '_theta_' Vr.tmmode],Tint);
stepTable = mms.db_get_ ts(datasetName,[pref '_steptable_parity_' Vr.tmmode],Tint);
res = irf.ts_skymap(dist.time,dist.data,[],phi.data,theta.data,'energy0',energy0.data,'energy1',energy1.data,'esteptable',stepTable.data);            

%% Fast data
tint = irf.tint('2015-10-16T12:01:00.00Z/2015-10-16T13:59:00.00Z');
Tint = tint;
c_eval('datasetName = ''mms?_fpi_fast_l2_dis-dist'';',ic)
c_eval('pref = ''mms?_dis'';',ic)
Vr.tmmode = 'fast';
dist = mms.db_get_ts(datasetName,[pref '_dist_' Vr.tmmode],Tint);
energy = mms.db_get_variable(datasetName,[pref '_energy_' Vr.tmmode],Tint(1));
phi = mms.db_get_variable(datasetName,[pref '_phi_' Vr.tmmode],Tint);
theta = mms.db_get_variable(datasetName,[pref '_theta_' Vr.tmmode],Tint);
res = irf.ts_skymap(dist.time,dist.data,[],phi.data,theta.data,'energy0',energy.data,'energy1',energy.data,'esteptable',TSeries(dist.time,zeros(dist.length,1)));            