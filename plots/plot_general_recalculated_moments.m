ic = 1;
tint = irf.tint('2017-08-20T01:53:00.00Z/2017-08-20T01:57:00.00Z'); %20160101090314
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');   



c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)

c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?); toc;',ic) 

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('tic; ne? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',tint); toc;',ic);
c_eval('tic; ni? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_dbcs_brst'',tint); toc;',ic);
c_eval('tic; gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?); toc;',ic)


c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)
%%
time = irf_time('2017-08-20T01:55:27.00Z','utc>EpochTT');

elim = [500 Inf];
tic; emoms_e500 = mms.psd_moments(ePDist1.elim(elim),scPot1); toc

elim = [100 Inf];
tic; emoms_e100 = mms.psd_moments(ePDist1.elim(elim),scPot1); toc

momVe1 = emoms_e100.V_psd;
momTe1 = emoms_e100.T_psd;
c_eval('momfacTe? = mms.rotate_tensor(momTe?,''fac'',gseB?); momfacTe?.units = ''eV''; momfacTe?.coordinateSystem = ''FAC'';',ic)

%% Figure: overview

fontsize = 12;

h = irf_plot(6);


if 1 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',fontsize);
end 


if 1 % Vi
  hca = irf_panel('Vi xyz');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
end

if 1 % Vi
  hca = irf_panel('Vi xyz parperp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x.tlim(tint),gseVi?perp.y.tlim(tint),gseVi?perp.z.tlim(tint),gseVi?par.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','||'},[0.98,0.3],'fontsize',fontsize);
end

if 0 % Ve
  hca = irf_panel('Ve xyz');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint)*1e-3,gseVe?.y.tlim(tint)*1e-3,gseVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
end

if 1 % Ve
  hca = irf_panel('Ve xyz partial');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{momVe?.x.tlim(tint)*1e-3,momVe?.y.tlim(tint)*1e-3,momVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
end
if 1 % Te
  hca = irf_panel('Te xyz partial');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{momfacTe?.xx,momfacTe?.yy,momfacTe?.zz},''comp'');',ic)  
  
  hca.YLabel.String = {'T_e (eV)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'||','\perp1','\perp2'},[0.98,0.3],'fontsize',fontsize);
end

if 1 % dEFlux
  hca = irf_panel('dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,ePDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  
  %hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  colormap(irf_colormap('thermal'))
  hca.YScale = 'log';
  hca.NextPlot = "add";
  plot(hca,hca.XLim,elim(1)*[1 1],'k--')
  hca.NextPlot = "replace";
end

h(end).XTickLabelRotation = 0;
irf_pl_mark(h,time,'black')
irf_zoom(h,'x',gseB1.time)
irf_plot_axis_align

%% Figure: overview, with distribution

time = irf_time('2017-08-20T01:55:26.90Z','utc>EpochTT');
fontsize = 12;

%h = irf_plot(5);
[h,h2] = initialize_combined_plot('leftright',5,1,1,0.7,'vertical');


if 1 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',fontsize);
end 


if 1 % Vi
  hca = irf_panel('Vi xyz');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
end

if 0 % Ve
  hca = irf_panel('Ve xyz');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint)*1e-3,gseVe?.y.tlim(tint)*1e-3,gseVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
end

if 1 % Ve
  hca = irf_panel('Ve xyz partial');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{momVe?.x.tlim(tint)*1e-3,momVe?.y.tlim(tint)*1e-3,momVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
end
if 1 % Te
  hca = irf_panel('Te xyz partial');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{momfacTe?.xx,momfacTe?.yy,momfacTe?.zz},''comp'');',ic)  
  
  hca.YLabel.String = {'T_e (eV)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'||','\perp1','\perp2'},[0.98,0.3],'fontsize',fontsize);
end

if 1 % dEFlux
  hca = irf_panel('dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,ePDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  
  %hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  colormap(irf_colormap('thermal'))
  hca.YScale = 'log';
  hca.NextPlot = "add";
  plot(hca,hca.XLim,elim(1)*[1 1],'k--')
  hca.NextPlot = "replace";
end

h(end).XTickLabelRotation = 0;
%irf_zoom(h,'x',gseB1.time)
irf_zoom(h,'x',time + 10*[-1 1])
irf_plot_axis_align
% irf_pl_mark(h,time,'black')

%%

for dt = -100:3:100
  time_tmp = time + dt*0.03;
  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h,time_tmp,'black');

  pdist = ePDist1.elim([100 Inf]).tlim(time_tmp+3*0.03*0.5*[-1 1]);
  vpar = mean(gseB1.tlim(time_tmp+3*0.03*0.5*[-1 1]).data,1);
  vpar = vpar/norm(vpar);
  perp1 = cross(vpar,cross([0 0 1],vpar));
  perp1 = perp1/norm(perp1);
  
  isub = 1;
  hca = h2(isub); isub = isub + 1;
  
  vdf = pdist.reduce('2D',perp1,vpar);
  vdf.plot_plane(hca,'contour',10)
  axis(hca,'square')
  hca.XLabel.String = 'v_{\perp1} (km/s)';
  hca.YLabel.String = 'v_{||} (km/s)';
  hca.XLim = 70*[-1 1];
  hca.YLim = 70*[-1 1];
  hca.XTick = hca.YTick;
  time_str = sprintf('%sT%s-%s',pdist.time(1).utc('yyyy-mm-dd'),pdist.time(1).utc('HH:MM:SS.mmm'),pdist.time(end).utc('HH:MM:SS.mmm'));
  hca.Title.String = time_str;
  cn.print(['psbl_', time_str])
end





