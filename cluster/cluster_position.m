downloadDirectory = '/Users/cno062/Data/Cluster';

tint = '2000-08-01T00:00:00.000Z/2001-12-31T24:00:00.000Z';
tint = '2003-01-01T05:00:00.000Z/2005-01-01T05:00:30.000Z';
%caa_download(tint,'C1_CP_FGM_FULL')
caa_download(tint,'C1_CP_FGM_5VPS',['downloadDirectory=',downloadDirectory])


%c_eval('gsePos?=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_SPIN'',''mat'');',3:4);         
%%

dobj = dataobj('./C1_CP_FGM_SPIN/C1_CP_FGM_SPIN__20000801_000000_20011201_000000_V170222.cdf');
R1_orig = irf.ts_vec_xyz(irf_time(dobj.data.time_tags__C1_CP_FGM_SPIN.data,'epoch>epochTT'),dobj.data.sc_pos_xyz_gse__C1_CP_FGM_SPIN.data);
units = irf_units;

%% Write data to text file
time = R1.time.utc;
x = R1.x.data;
y = R1.y.data;
z = R1.z.data;
datatable = table(time,x,y,z);
writetable(datatable,'Cluster_position_2000-12-03_2021_11_30.txt');
%%
R1 = R1_orig.resample(R1_orig.time(1):60:R1_orig.time(end));
R1 = R1_orig;
fontsize = 20;
nrows = 3;
ncols = 1;
%h = setup_subplots(nrows,ncols);
h(1) = subplot(4,1,1);
h(2) = subplot(4,1,2:3);
h(3) = subplot(4,1,4);
isub = 1;

hca = h(isub); isub = isub + 1;
irf_plot(hca,R1*1e3/units.RE,'linewidth',1)
hca.YLabel.String = 'Position (R_E)';
hca.YLabel.Interpreter = 'tex';
irf_legend(hca,{'x','y','z'},[0.02 0.1],'fontsize',fontsize)
irf_zoom(hca,'x',R1.time)
%
hca = h(isub); isub = isub + 1;
scatter(hca,R1.x.data*1e3/units.RE,R1.y.data*1e3/units.RE,3,(R1.time-R1.time(1))/60/60/24)
hca.XLabel.String = 'X_{GSE} (R_E)';
hca.YLabel.String = 'Y_{GSE} (R_E)';
hca.XLabel.Interpreter = 'tex';
hca.YLabel.Interpreter = 'tex';
hcb = colorbar('peer',hca);
hcb.YLabel.String = ['Days since ' irf_time(R1.time(1).utc,'utc>utc_yyyy-mm-dd')];
axis(hca,'square')
axis(hca,'equal')
hca.XGrid = 'on';
hca.YGrid = 'on';

hca = h(isub); isub = isub + 1;
scatter(hca,R1.x.data*1e3/units.RE,R1.z.data*1e3/units.RE,3,(R1.time-R1.time(1))/60/60/24)
hca.XLabel.String = 'X_{GSE} (R_E)';
hca.YLabel.String = 'Z_{GSE} (R_E)';
hca.XLabel.Interpreter = 'tex';
hca.YLabel.Interpreter = 'tex';
hcb = colorbar('peer',hca);
hcb.YLabel.String = ['Days since ' irf_time(R1.time(1).utc,'utc>utc_yyyy-mm-dd')];
axis(hca,'equal')
hca.XGrid = 'on';
hca.YGrid = 'on';

h(2).XLim = h(3).XLim;
c_eval('h(?).FontSize = fontsize;',1:numel(h));
