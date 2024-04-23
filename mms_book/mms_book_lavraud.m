localuser = 'cno062';
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db'); 

tint = irf.tint('2015-10-16T10:32:30.00Z/2015-10-16T10:34:10.00Z');
ic = 4;
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

%%
tint_plot = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:38.00Z');
elim = [165 210];

ePitch4 = ePDist4(1:1:ePDist4.length).tlim(tint_plot + [-1 1]).pitchangles(dmpaB4,0:7.5:180);
specrec = ePitch4.elim(elim).tlim(tint_plot + [-1 1]).specrec;
specrec.p = smooth2(specrec.p,1);
tint_angle = irf.tint('2015-10-16T10:33:25.00Z/2015-10-16T10:33:38.00Z');
tsAngle = trapping_angle(30,90,gseB4.abs);
tsAngle2 = trapping_angle(20,90,gseB4.abs);

%h = irf_plot(2);
[h,h2] = initialize_combined_plot(2,2,2,0.55,'vertical');

hca = irf_panel('Babs');
set(hca,'colororder',mms_colors('1'))
irf_plot(hca,dmpaB4.abs,'linewidth',1)
hca.YLabel.String = '|B| (nT)';
hca.XGrid = 'off';
hca.YGrid = 'off';


hca = irf_panel('pitchangles');
set(hca,'colororder',mms_colors('1'))
[hsp,hcb] = irf_spectrogram(hca,specrec,'log');
hcb.YLabel.String = 'f_e (s^3/km^6)';
colormap(hca,pic_colors('waterfall'))
colormap(hca,pic_colors('candy4'))
%colormap(hca,'jet')
hold(hca,'on')
irf_plot(hca,tsAngle.tlim(tint_angle),'k--','linewidth',1)
%irf_plot(hca,tsAngle2.tlim(tint_angle+[1.2 -1.7]),'k:','linewidth',1)
hold(hca,'off')
hca.XTickLabelRotation = 0;
hca.CLim = [-28.2 -27];
hca.Box = 'on';
hca.Layer = 'top';
hca.XGrid = 'off';
hca.YGrid = 'off';
hca.YLabel.String = '\theta (deg)';
hca.YLabel.Interpreter = 'tex';

irf_zoom(h,'x',tint_plot)
irf_plot_axis_align
c_eval('h(?).FontSize = 16;',1:numel(h))

% 2D reduced plots
vg = -20000:500:20000;
vlim = 0.99*15*[-1 1];

isub = 1;

if 1 % 1
  hca = h2(isub); isub = isub + 1;
  time = irf_time('2015-10-16T10:33:24.50Z','utc>EpochTT');
  time_int = time+0.5*0.03*[-1 1];
  scpot = scPot4.resample(ePDist4);
  par = mean(dmpaB4.tlim(time_int).data,1); par = par/norm(par);
  perp1 = mean(dslE4.tlim(time_int).data,1); perp1 = perp1/norm(perp1);
  perp1 = cross(par,cross(perp1,par));
  perp2 = cross(par,perp1);
  pdist = ePDist4.tlim(time_int);
  vdf = pdist.reduce('2D',par,perp1,'vg',vg,'scpot',scpot);
  vdf.plot_plane(hca,'contour',5);
  colormap(hca,pic_colors('candy4'))
  hca.XLim = vlim;
  hca.YLim = vlim;
  irf_pl_mark(h,time,'k','linestyle',':')  
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{BxExB} (10^3 km/s)';
end 
if 1 % 3
  hca = h2(isub); isub = isub + 1;
  time = irf_time('2015-10-16T10:33:27.50Z','utc>EpochTT');
  time_int = time+0.5*0.03*[-1 1];
  scpot = scPot4.resample(ePDist4);
  par = mean(dmpaB4.tlim(time_int).data,1); par = par/norm(par);
  perp1 = mean(dslE4.tlim(time_int).data,1); perp1 = perp1/norm(perp1);
  perp1 = cross(par,cross(perp1,par));
  perp2 = cross(par,perp1);
  pdist = ePDist4.tlim(time_int);
  vdf = pdist.reduce('2D',par,perp1,'vg',vg,'scpot',scpot);
  vdf.plot_plane(hca,'contour',5);
  colormap(hca,pic_colors('candy4'))
  hca.XLim = vlim;
  hca.YLim = vlim;
  irf_pl_mark(h,time,'k','linestyle',':')
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{BxExB} (10^3 km/s)';
end
if 1 % 2
  hca = h2(isub); isub = isub + 1;
  time = irf_time('2015-10-16T10:33:26.00Z','utc>EpochTT');
  time_int = time+0.5*0.03*[-1 1];
  scpot = scPot4.resample(ePDist4);
  par = mean(dmpaB4.tlim(time_int).data,1); par = par/norm(par);
  perp1 = mean(dslE4.tlim(time_int).data,1); perp1 = perp1/norm(perp1);
  perp1 = cross(par,cross(perp1,par));
  perp2 = cross(par,perp1);
  pdist = ePDist4.tlim(time_int);
  vdf = pdist.reduce('2D',par,perp1,'vg',vg,'scpot',scpot);
  vdf.plot_plane(hca,'contour',5);
  colormap(hca,pic_colors('candy4'))
  hca.XLim = vlim;
  hca.YLim = vlim;
  irf_pl_mark(h,time,'k','linestyle',':')
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{BxExB} (10^3 km/s)';
end
if 1 % 4
  hca = h2(isub); isub = isub + 1;
  time = irf_time('2015-10-16T10:33:29.70Z','utc>EpochTT');
  %time = irf_time('2015-10-16T10:33:26.443Z','utc>EpochTT');
  time_int = time+2*0.5*0.03*[-1 1];
  scpot = scPot4.resample(ePDist4);
  par = mean(dmpaB4.tlim(time_int).data,1); par = par/norm(par);
  perp1 = mean(dslE4.tlim(time_int).data,1); perp1 = perp1/norm(perp1);
  perp1 = cross(par,cross(perp1,par));
  perp2 = cross(par,perp1);
  pdist = ePDist4.tlim(time_int);
  vdf = pdist.reduce('2D',par,perp1,'vg',vg,'scpot',scpot);
  vdf.plot_plane(hca,'contour',5);
  colormap(hca,pic_colors('candy4'))
  hca.XLim = vlim;
  hca.YLim = vlim;
  irf_pl_mark(h,time,'k','linestyle',':')
  hca.XLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.String = 'v_{BxExB} (10^3 km/s)';
end
if 0
  hca = h2(isub); isub = isub + 1;
  time = irf_time('2015-10-16T10:33:26.60Z','utc>EpochTT');
  time_int = time+0.5*0.03*[-1 1];
  scpot = scPot4.resample(ePDist4);
  par = mean(dmpaB4.tlim(time_int).data,1); par = par/norm(par);
  perp1 = mean(dslE4.tlim(time_int).data,1); perp1 = perp1/norm(perp1);
  perp1 = cross(par,cross(perp1,par));
  perp2 = cross(par,perp1);
  pdist = ePDist4.tlim(time_int);
  vdf = pdist.reduce('2D',par,perp1,'vg',vg,'scpot',scpot);
  vdf.plot_plane(hca,'contour',2);
  colormap(hca,pic_colors('candy4'))
  hca.XLim = vlim;
  hca.YLim = vlim;
  irf_pl_mark(h,time,'k','linestyle',':')
end
if 0
  hca = h2(isub); isub = isub + 1;
  time = irf_time('2015-10-16T10:33:30.260Z','utc>EpochTT');
  time_int = time+0.5*0.03*[-1 1];
  scpot = scPot4.resample(ePDist4);
  par = mean(dmpaB4.tlim(time_int).data,1); par = par/norm(par);
  perp1 = mean(dslE4.tlim(time_int).data,1); perp1 = perp1/norm(perp1);
  perp1 = cross(par,cross(perp1,par));
  perp2 = cross(par,perp1);
  pdist = ePDist4.tlim(time_int);
  vdf = pdist.reduce('2D',perp2,perp1,'vg',vg,'scpot',scpot,'vint',[-2000 2000]);
  vdf.plot_plane(hca);
  colormap(hca,pic_colors('candy4'))
  hca.XLim = vlim;
  hca.YLim = vlim;
  irf_pl_mark(h,time,'k','linestyle',':')
end

compact_panels(h2,0.02,0.01)

linkprop(h2,{'XLim','YLim','CLim'})
h2(1).CLim = [-12 -6.4];
c_eval('axis(h2(?),''square'')',1:4)
