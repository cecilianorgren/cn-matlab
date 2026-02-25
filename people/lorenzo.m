%% Define time interval
tint = irf.tint('2015-10-16T10:30:00.00Z/2015-10-16T10:30:40.00Z');



dobj = dataobj([db_info.local_file_db_root '/mms1/fpi/brst/l2/des-dist/2015/10/27/' 'mms1_fpi_brst_l2_des-dist_20151027123344_v2.1.0.cdf']);
desDist1 = mms.variable2ts(get_variable(dobj,'mms1_des_dist_brst'));

time = irf_time('2015-10-16T10:30:00.00Z','utc>epochtt');
mms.plot_projection(hca,desDist1,'tint',time);
mms.plot_projection(hca,desDist1,'tint',time,'xyz',[1 0 0;0 0 1; 0 -1 0]);

mms.plot_skymap(desDist1,'tint',time,'energy',100)

%%

figure
nRows = 4;
nCols = 1;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end

isub = 1;

hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',time,'energy',100)

hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,desDist1,'tint',time,'energy',200)

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,desDist1,'tint',time,'xyz',[1 0 0;0 1 0; 0 0 1]);

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,desDist1,'tint',time,'xyz',[1 0 0;0 0 1; 0 -1 0]);

%%
c_eval('B0 = gseB?.resample(time).data;',ic); 
z = double(irf_norm(B0));
y = cross(z,[1 0 0]);
x = cross(y,z)
%%
gseB1=mms.db_get_ts('mms1_fgm_brst_l2','mms1_fgm_b_gse_brst_l2',tint);
gseB2=mms.db_get_ts('mms2_fgm_brst_l2','mms2_fgm_b_gse_brst_l2',tint);
gseB3=mms.db_get_ts('mms3_fgm_brst_l2','mms3_fgm_b_gse_brst_l2',tint);
gseB4=mms.db_get_ts('mms4_fgm_brst_l2','mms4_fgm_b_gse_brst_l2',tint);

% OR

c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); gseB? = gseB?.clone(gseB?.time,gseB?.data(:,1:3));',1:4);