mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
tint = irf.tint('2015-10-16T10:32:40.00Z/2015-10-16T10:34:10.00Z');
ic = 1;
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
%c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('tic; dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?); toc;',ic)

c_eval('gseVExB? = cross(dslE?,dmpaB?.resample(dslE?.time))/dmpaB?.abs.resample(dslE?.time)/dmpaB?.abs.resample(dslE?.time)*1e3;',ic) % km/s
%%

time = irf_time('2015-10-16T10:33:30.326Z','utc>EpochTT');
tint_dist = time+[-0.015 0.015] + 0*0.03;
dist = ePDist1.tlim(tint_dist);
%scpot = mean(scPot1.tlim(tint_dist).data,1);
scpot = scPot1.resample(dist)

par = mean(dmpaB1.tlim(tint_dist).data,1); % B
par = par/norm(par);
perp1 = mean(dslE1.tlim(tint_dist).data,1); % E
perp1 = perp1/norm(perp1);
perp2 = cross(par,perp1); % ExB

vint = inf*[-4000 4000];
vg = -10000:100:10000;

f2D = dist.reduce('2D',perp1,perp2,'vg',vg,'nMC',800,'scpot',scpot,'lowerelim',10,'vint',vint);
[h_surf,h_axis,h_all] = f2D.plot_plane('contour',[1 1]*10^(-8.2)); % 10.^(-10:0.5:-5)
colormap(pic_colors('candy4'))
h_axis.YLim = [min(vg) max(vg)]*1e-3;
h_axis.XLim = [min(vg) max(vg)]*1e-3;
%h_axis.CLim = [-9 -6];
h_axis.CLim = [-9.8 -6.5];
h_axis.XGrid = 'on';
h_axis.YGrid = 'on';
h_axis.Layer = 'top';
h_axis.XTick = -10:2:10;
h_axis.YTick = -10:2:10;
h_axis.XLabel.String = 'v_{E} (10^3 km/s)';
h_axis.YLabel.String = 'v_{ExB} (10^3 km/s)';
h_axis.FontSize = 16;
disp('done')
%[h_surf,h_axis,h_all] = f2D(100).plot_plane;