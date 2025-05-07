%% Set up database
units = irf_units;
irf.log('critical')
ic_dist = 1;
ic = 1;

localuser = datastore('local','user');
localuser = 'cecilia';
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
%mms.db_init('local_file_db',['/Volumes/DataRaid/MMS']);
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db'); 

%% Time intervals
%tint = irf.tint('2017-07-06T15:45:45.00Z/2017-07-06T22:48:45.00Z');
tint = irf.tint('2017-07-03T05:26:13.00Z/2017-07-03T05:27:03.00Z'); % vertical crossing, come cold ions I think
tint = irf.tint('2017-07-03T21:54:03.00Z/2017-07-03T21:55:53.00Z'); % only skimming, cold ions
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z'); % cold ions, my 2020 separatrix paper
tint = irf.tint('2017-07-06T08:16:03.00Z/2017-07-06T08:18:13.00Z'); % cold ions, skimming
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z'); % cold ions, mostly skimming, my EDI paper
tint = irf.tint('2017-07-06T14:06:53.00Z/2017-07-06T14:08:53.00Z'); % cold ions, skimming
tint = irf.tint('2017-07-06T14:08:53.00Z/2017-07-06T14:10:20.00Z'); % cold ions, skimming
tint = irf.tint('2017-07-06T15:28:03.00Z/2017-07-06T15:30:03.00Z'); % cold ions, a few vertical crossings
tint = irf.tint('2017-07-06T15:41:33.00Z/2017-07-06T15:44:03.00Z'); % cold ions, some almost crossings
tint = irf.tint('2017-07-06T15:46:23.00Z/2017-07-06T15:48:53.00Z'); % cols ions, some almost crossings
tint = irf.tint('2017-07-06T16:02:13.00Z/2017-07-06T16:04:23.00Z'); % cold ions, some full crossings
tint = irf.tint('2017-07-06T16:35:43.00Z/2017-07-06T16:38:13.00Z'); % cold ions, some full and partial crossings
tint = irf.tint('2017-07-17T15:04:43.00Z/2017-07-17T15:05:23.00Z'); % cold ions, low B overall, hard to interepret perhaps
tint = irf.tint('2017-07-18T01:36:43.00Z/2017-07-18T01:38:23.00Z'); % cold ions, skimming, low fluxes?


tint = irf.tint('2017-07-06T15:28:03.00Z/2017-07-06T15:30:03.00Z'); % cold ions, a few vertical crossings

tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); % torbert science event


str_tint = tint(1).utc('yyyymmdd_HHMMSS');
%% Load data
disp('Loading B, E, Vsc...')
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
%c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',1:4);
% Electric field
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
%c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',1:4);
%c_eval('gsmE? = c_coord_trans(''GSE'',''GSM'',gseE?);',1:4)

% Spacecraft potential
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

disp('Loading n, v...')
% Density
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
%c_eval('gsmVe? = c_coord_trans(''GSE'',''GSM'',gseVe?);',ic)
%c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?);',ic)

disp('Loading eVDFs...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic_dist) % missing some ancillary data
disp('Loading iVDFs...')
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic_dist) % missing some ancillary data
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic_dist) % missing some ancillary data

disp('Calculating derived quantities...')
% Remove all one-count "noise"

%c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic_dist)
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?_counts.data./iPDistErr?.data).^2;',ic_dist)
%c_eval('iPDist?.data(iPDist?.data < iPDistErr?.data*1.01) = 0;',ic)


c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('gseVExB?_EyBz = gseE?.y.resample(gseB?.time).*gseB?.z/gseB?.abs/gseB?.abs*1e3; gseVExB?_EyBz.units = '''';',ic) % km/s
c_eval('gseVExB?_EzBy = -1*gseE?.z.resample(gseB?.time).*gseB?.y/gseB?.abs/gseB?.abs*1e3; gseVExB?_EzBy.units = '''';',ic) % km/s
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('[gseVi?par,gseVi?perp] = irf_dec_parperp(gseB?,gseVi?); gseVi?par.name = ''Vi par''; gseVi?perp.name = ''Vi perp'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.name = ''v_ixB'';',ic)
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.name = ''v_exB'';',ic)


%% Make reduced VDFs
ic = ic_dist;
disp('Preparing reduced distributions.')
vint = [-Inf Inf];
iELim = [200 40000]; % for tail



%c_eval('if1Dx?_700 = iPDist?.elim([700 Inf]).reduce(''1D'',[1 0 0],''vint'',vint);',ic)
%c_eval('if1Dy?_700 = iPDist?.elim([700 Inf]).reduce(''1D'',[0 1 0],''vint'',vint);',ic)
%c_eval('if1Dz?_700 = iPDist?.elim([700 Inf]).reduce(''1D'',[0 0 1],''vint'',vint);',ic)

c_eval('if1Dx? = iPDist?.elim(iELim).reduce(''1D'',[1 0 0],''vint'',vint,''scpot'',scPot?);',ic)
c_eval('if1Dy? = iPDist?.elim(iELim).reduce(''1D'',[0 1 0],''vint'',vint,''scpot'',scPot?);',ic)
c_eval('if1Dz? = iPDist?.elim(iELim).reduce(''1D'',[0 0 1],''vint'',vint,''scpot'',scPot?);',ic)


%% Figure: Overview 1
h = irf_plot(9);

nMovMean = 5;
fontsize = 12;

if 1 % B GSE
  hca = irf_panel('B GSE');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',fontsize);
end 
if 0 % B LMN
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
end 
if 1 % E GSE
  hca = irf_panel('E GSE');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('toplot = gseE?.resample(gseVi?);',ic)
  %c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  c_eval('irf_plot(hca,{toplot.x,toplot.y,toplot.z},''comp'');',ic)
  hca.YLabel.String = {'E (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',fontsize);
end 

if 1 % Ve, GSE
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.98],'fontsize',fontsize);
end
if 0 % Ve, LMN
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % Vi, GSE
  hca = irf_panel('Vi GSE');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % Vi, LMN
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % Pi
  hca = irf_panel('Pi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPi?.xx.tlim(tint),mvaPi?.yy.tlim(tint),mvaPi?.zz.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'P_i (nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % Pi
  hca = irf_panel('Pi off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPi?.xy.tlim(tint),mvaPi?.xz.tlim(tint),mvaPi?.yz.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'P_i (nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % VExB, GSE
  hca = irf_panel('ExB GSE');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.resample(gseVi?),gseVExB?_EyBz.resample(gseVi?),gseVExB?_EzBy.resample(gseVi?)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_{ExB,x} (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{ExB,x}','E_yB_z','-E_zB_y'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % dEFlux ion
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  hca.NextPlot = "add";
  c_eval('irf_plot(hca,irf.ts_scalar(tint,iELim(1)*[1 1]),''k--'')',ic) 
  hca.NextPlot = "replace";
end
if 0 % dEFlux ion, movmean, remove onecounts
  hca = irf_panel('ion dEF omni movmean rem onecounts');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.1],'fontsize',fontsize,'color','k');
end

if 1 % fi red x
  hca = irf_panel('fi x');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = if1Dx?.specrec;',ic)
  %specrec.f = specrec.f - v_xline;
  c_eval('irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
    
  hca.NextPlot = "add";
  c_eval('irf_plot(hca,gseVExB?.resample(gseVi?).x,''k-'')',ic)
  %c_eval('irf_plot(hca,irf.ts_scalar(mvaVi?.time,repmat(v_xline,[mvaVi?.length,1])),''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{ix} (km/s)'};

  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 1 % fi red y
  hca = irf_panel('fi y');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,if1Dy?.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  c_eval('irf_plot(hca,gseVExB?.resample(gseVi?).y,''k-'')',ic)
  %c_eval('irf_plot(hca,mvaVi?.y,''k-'')',ic)
  hca.NextPlot = "replace";
  
  hca.YLabel.String = {'v_{iy} (km/s)'};
end
if 1 % fi red z
  hca = irf_panel('fi z');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,if1Dz?.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  c_eval('irf_plot(hca,gseVExB?.resample(gseVi?).z,''k-'')',ic)
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iz} (km/s)'};
end

if 0 % fi red L
  hca = irf_panel('fi L');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_L.specrec;',ic)
  %specrec.f = specrec.f - v_xline;
  c_eval('irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
    
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.x,''k-'')',ic)
  %c_eval('irf_plot(hca,irf.ts_scalar(mvaVi?.time,repmat(v_xline,[mvaVi?.length,1])),''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iL} (km/s)'};

  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 0 % fi red M
  hca = irf_panel('fi M');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_M.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.y,''k-'')',ic)
  hca.NextPlot = "replace";
  
  hca.YLabel.String = {'v_{iM} (km/s)'};
end
if 0 % fi red L
  hca = irf_panel('fi N');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_N.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iN} (km/s)'};
end
if 0 % max shear direction
  hca = irf_panel('fi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_e1.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iN} (km/s)'};
end
if 0 % e1
  hca = irf_panel('Vi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{e1.x,e1.y,e1.z},''comp'');',ic)  
  
  hca.YLabel.String = {'e_1 (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end


%irf_zoom(h,'x',tint_figure)
irf_zoom(h,'x',tint)
%irf_zoom(h(1:3),'y')
%irf_pl_mark(h,time_xline,'black','linestyle',':')
%irf_pl_mark(h,time_xline_ion,'red','linestyle',':')
%colormap(irf_colormap('thermal'))
colormap([pic_colors('candy_gray'); 1 1 1])
colormap([pic_colors('candy_gray')])
colormap([pic_colors('candy4')])
irf_plot_axis_align
h(end).XTickLabelRotation = 0;
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
c_eval('h(?).Layer = ''top''; h(?).Box = ''on'';',1:numel(h))
%hlinks1 = linkprop(h(3:4),{'CLim'});
%hlinks2 = linkprop(h(5:7),{'CLim'});
%hlinks2 = linkprop(h(4:6),{'CLim'});
%h(1).Title.String = sprintf('N_{movmean} = %g',nMovMean);
c_eval('h(?).FontSize = 14;',1:numel(h))


if 1
  %%
  linkprop(h((numel(h)-2):end),{'YLim'})
  linkprop(h((numel(h)-2):end),{'CLim'})
  h(end).CLim = [-4 -1];
  h(end).YLim = [-600 600];
  h(end).YLim = [-800 800];
  h(end).YLim = 2*[-1000 1000];
end

%% Figure: 3D isosurface with 3D pick-up surface

units = irf_units;
ic = 2;
times_utc = ['2015-10-16T13:07:02.160Z';...
             '2015-10-16T13:07:02.190Z';...
             '2015-10-16T13:07:02.220Z';...
             '2015-10-16T13:07:02.250Z';...
             '2015-10-16T13:07:02.280Z'];
times = irf_time(times_utc,'utc>EpochTT');

%times = EpochTT(['2017-07-11T22:34:01.300000000Z';
%  '2017-07-11T22:34:01.800000000Z';...
%  '2017-07-11T22:34:02.300000000Z';...
%  '2017-07-11T22:34:02.800000000Z';...
%  '2017-07-11T22:34:03.300000000Z']);
%times = times + 0.090;
times = times(1:3);
nt = times.length;

vint_red = [-1000 1000];
vint = [3500 5000];
vint = [3000 5000] + 500;
eint = [50 90];
eint = [50 90];
vlim_red = [-10 10];
eint = units.me*(vint*1e3).^2/2/units.eV;
c_eval('PDist = ePDist?;',ic)
c_eval('dmpaB = dmpaB?;',ic)
c_eval('dslE = dslE?;',ic)
c_eval('gseE = gseE?;',ic)
c_eval('scPot = scPot?;',ic)


%if exist('h','var'); delete(h)

nrows = 1;
ncols = nt;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;
%nt = 1;

set(gcf,'defaultLegendAutoUpdate','off');

for it = 1:nt % 3D isosurface
    
  hca = h(isub); isub = isub + 1;
  time = times(it);
  elimit = 40;
  elimit = 35;
  elimit = 30;
  pdist = PDist.elim([elimit Inf]).tlim(time+0.499*0.03*[-1 1]);  
  pdist = PDist.tlim(time+0.499*0.03*[-1 1]);  
  pdist.data(:,1:5,:,:) = 3e-27;


  scP = scPot.tlim(time+0.5*0.03*[-1 1]); scP = mean(scP.data,1);
  B = dmpaB.tlim(time+0.5*0.03*[-1 1]); B = mean(B.data,1); Bn = B/norm(B);
  E = dslE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); En = E/norm(E);
  %E = gseE.tlim(time+0.5*0.03*[-1 1]); E = mean(E.data,1); En = E/norm(E);
  Eperp = cross(Bn,cross(E,Bn)); Eperpn = Eperp\norm(Eperp);
  ExB = cross(E,B); ExBn = Eperp\norm(Eperp);
  
  nSmooth = 3;
  
  iso_values = [0.5e-27, 6e-27];
  %iso_values = [1.5e-30 22e-30];
  iso_values = iso_values(2); 
  set(hca,'colororder',[1 0.5 0.5; 0.5 1 0.5]) 
  set(hca,'colororder',[0.2 0.5 1; 1 0.5 0.5; 0.5 1 0.5]) 
  %Tgse = [Lgse;Mgse;Ngse];
  c_eval('Tdsl = [tsLdsl?.resample(pdist).data;tsMdsl?.resample(pdist).data;tsNdsl?.resample(pdist).data];',ic)
  hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',Tdsl);
  %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',Tdsl);
  
  hs.Patch(1).FaceAlpha = 0.9;
  %hs.Patch(1).FaceAlpha = 0.3;
  %hs.Patch(2).FaceAlpha = 0.9;
  
  legs = arrayfun(@(x)sprintf('%g',x),hs.Values,'UniformOutput',false);
  hleg = legend(hs.Patch,legs,'location','northeast');
  hleg.Title.String = sprintf('f_%s (%s)',pdist.species(1),pdist.units);
  
  
  if 1 % plot E, ExB, B    
    hold(hca,'on')
    qmax = 10000; 
    nE = E/norm(E);
    nEperp = Eperp/norm(Eperp);
    nB = B/norm(B);
    nExB = ExB/norm(ExB);
    colors_quiv_B = [0 0 0];
    colors_quiv_E = [1 0 0];
    colors_quiv_Eperp = [1 0.5 0.5];
    colors_quiv_ExB = [0 0.5 1]; 

    if 1
      nE = nE*hs.Rotation';
      nExB = nExB*hs.Rotation';
      nB = nB*hs.Rotation';
      nEperp = nEperp*hs.Rotation';
    end
    
    fontsize = 13;
    fontweight = 'bold';
    quiver3(hca,qmax*nB(1)*-1,qmax*nB(2)*-1,qmax*nB(3)*-1,qmax*nB(1)*2,qmax*nB(2)*2,qmax*nB(3)*2,'color',colors_quiv_B,'linewidth',2)
    %text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,'B','color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nB(1)*1,qmax*nB(2)*1,qmax*nB(3)*1,sprintf('B=%.0fnT',norm(B)),'color',colors_quiv_B,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*nE(1)*-1,qmax*nE(2)*-1,qmax*nE(3)*-1,qmax*nE(1)*2,qmax*nE(2)*2,qmax*nE(3)*2,'color',colors_quiv_E,'linewidth',2)
    %text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,'E','color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nE(1)*1,qmax*nE(2)*1,qmax*nE(3)*1,sprintf('E=%.0fmV/m',norm(E)),'color',colors_quiv_E,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nEperp(1)*-1,qmax*nEperp(2)*-1,qmax*nEperp(3)*-1,qmax*nEperp(1)*2,qmax*nEperp(2)*2,qmax*nEperp(3)*2,'color',colors_quiv_Eperp,'linewidth',2)
    %text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,'E_\perp','color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
    text(hca,qmax*nEperp(1)*1,qmax*nEperp(2)*1,qmax*nEperp(3)*1,['E_\perp' sprintf('=%.0fmV/m',norm(Eperp))],'color',colors_quiv_Eperp,'fontweight',fontweight,'fontsize',fontsize)
    
    quiver3(hca,qmax*nExB(1)*-1,qmax*nExB(2)*-1,qmax*nExB(3)*-1,qmax*nExB(1)*2,qmax*nExB(2)*2,qmax*nExB(3)*2,'color',colors_quiv_ExB,'linewidth',2)    
    text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,'E\times{B}','color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)
    %text(hca,qmax*nExB(1)*1,qmax*nExB(2)*1,qmax*nExB(3)*1,sprintf('E\times{B}=%.0fmV/m',norm(Eperp)),'color',colors_quiv_ExB,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')
    %set(hca,'colororder',colors_quiv)
    %irf_legend(hca,{'B','E','E_\perp','E\times{B}'},[0.98 0.98])
  end
  if 0 % plot LMN (obs in GSE!!)
    %%
     if 1
      L_ = L*hs.Rotation';
      M_ = M*hs.Rotation';
      N_ = N*hs.Rotation';
      
     end

    hold(hca,'on')
    colors_quiv_LMN = [0.5 0.5 0.5];
    quiver3(hca,qmax*L_(1)*-1,qmax*L_(2)*-1,qmax*L_(3)*-1,qmax*L_(1)*2,qmax*L_(2)*2,qmax*L_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*L_(1)*1,qmax*L_(2)*1,qmax*L_(3)*1,'L','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*M_(1)*-1,qmax*M_(2)*-1,qmax*M_(3)*-1,qmax*M_(1)*2,qmax*M_(2)*2,qmax*M_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*M_(1)*1,qmax*M_(2)*1,qmax*M_(3)*1,'M','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    quiver3(hca,qmax*N_(1)*-1,qmax*N_(2)*-1,qmax*N_(3)*-1,qmax*N_(1)*2,qmax*N_(2)*2,qmax*N_(3)*2,'color',colors_quiv_LMN,'linewidth',2)
    text(hca,qmax*N_(1)*1,qmax*N_(2)*1,qmax*N_(3)*1,'N','color',colors_quiv_LMN,'fontweight',fontweight,'fontsize',fontsize)

    hold(hca,'off')

  end
end

for ip = 1:numel(h)
  h(ip).XLabel.String = 'v_L (km/s)';
  h(ip).YLabel.String = 'v_M (km/s)';
  h(ip).ZLabel.String = 'v_N (km/s)';
  vlim = [-9 9]*1e3;
  vlim = 10*[-1 1]*1e3;
  h(ip).XLim = vlim;
  h(ip).YLim = vlim;
  h(ip).ZLim = vlim;
  h(ip).XTick = [-10:5:10]*1e3;
  h(ip).YTick = [-10:5:10]*1e3;
  h(ip).ZTick = [-10:5:10]*1e3;
  h(ip).XGrid = 'on';
  h(ip).YGrid = 'on';
  h(ip).ZGrid = 'on';
  h(ip).Box = 'on';
end

hlinks = linkprop(h,{'view'});
c_eval('h(?).FontSize = 12;',1:numel(h))

%view([0 0 1])

if 1 % Reset camlight
  %%
  hlight = findobj(h,'type','light');
  delete(hlight)
  for ip = 1:nt
    %camlight(h(ip),45,45);
    camlight(h(ip),45,10);
    camlight(h(ip),45,190);
  end
end

if 0 % Put axis at origin, missing labels, i.e. v_L (km/s)
  %%
  set(gca,'Visible','off')
  %c_eval('h(?).XTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).YTick = [-10:2:10]*1e3;',1:numel(h))
  %c_eval('h(?).ZTick = [-10:2:10]*1e3;',1:numel(h))
  hold(gca,'on')
  % DRAW AXIS LINEs
  plot3(get(gca,'XLim'),[0 0],[0 0],'k');
  plot3([0 0],[0 0],get(gca,'ZLim'),'k');
  plot3([0 0],get(gca,'YLim'),[0 0],'k');
  % GET TICKS
  X=get(gca,'Xtick');
  Y=get(gca,'Ytick');
  Z=get(gca,'Ztick');
  % GET LABELS
  XL=get(gca,'XtickLabel');
  YL=get(gca,'YtickLabel');
  ZL=get(gca,'ZtickLabel');
  % REMOVE TICKS
  set(gca,'Xtick',[]);
  set(gca,'Ytick',[]);
  set(gca,'Ztick',[]);
  % GET OFFSETS
  Xoff=diff(get(gca,'XLim'))./30;
  Yoff=diff(get(gca,'YLim'))./30;
  Zoff=diff(get(gca,'ZLim'))./30;
  % DRAW TICKS
  %%%%%%% THIS COULD BE VECTORiZeD %%%%%%%
  for i=1:length(X)
     plot3([X(i) X(i)],[0 0],[-Zoff Zoff],'k');
  end;
  for i=1:length(Y)
     plot3([-Xoff Xoff],[Y(i) Y(i)],[0 0],'k');
  end;
  for i=1:length(Z)
     plot3([-Xoff Xoff],[0 0],[Z(i) Z(i)],'k');
  end;
  % DRAW LABELS
  text(X,zeros(size(X)),zeros(size(X))-3.*Zoff,XL);
  text(zeros(size(Y))-3.*Xoff,Y,zeros(size(Y)),YL);
  text(zeros(size(Z))-3.*Xoff,zeros(size(Z)),Z,ZL);

  xlim_ = get(gca,'XLim');
  ylim_ = get(gca,'YLim');
  zlim_ = get(gca,'ZLim');
  text(gca,xlim_(2),0,0,'L (km/s)')
  text(gca,0,ylim_(2),0,'M (km/s)')
  text(gca,0,0,zlim_(2),'N (km/s)')
  
  hold(gca,'off')
end
%% 3D isosurface, to illustrate 3D tilt



elim = [2000 Inf];
nMovMean = 7;
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts);',ic)
pdist_all.data(:,1:20,:,:) = 0;
%pdist_all = pdist_all.elim(elim);

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
dt = -8;
time = time + dt;


c_eval('vExB = gseVExB?.tlim(time+0.150*0.5*nMovMean*[-1 1]).data;',ic)
vExB_mean = mean(vExB,1);
vExB_std = std(vExB,1);

h = setup_subplots(1,1);
isub = 1;

pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
%iso_values = [6.5e-27];
iso_values = [9e-27];
% adaptive
iso_values = prctile(pdist.data(:),[98.5]);
hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill');
%hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
%hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
axis(hca,'square')
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];
hca.ZLim = vlim*[-1 1];
%camlight(gca,90,-45)
%view(hca,[-1 -1 0.5])
%view(hca,[1 0.2 0.2])
view(hca,[0 1 0.2])
camlight(gca,0,0)
camlight(gca,180,0)
camlight(gca,0,90)
camlight(gca,0,-90)

if 1 % Add constant surface energy in ExB frame
  %
    vExB_mean(3) = 0;


  Ux_in_ExB_frame = 15000; % eV
  Uy_in_ExB_frame = 15000; % eV
  Uz_in_ExB_frame = 15000; % eV
  vabs_ExB = sqrt(sum(vExB_mean.^2));
  U0 = units.mp.*(vabs_ExB*1e3).^2/2/units.eV;
  Ux_in_ExB_frame = U0; % eV
  Uy_in_ExB_frame = U0; % eV
  Uz_in_ExB_frame = U0; % eV
  %Emax = 32000;
  vx_in_ExB = sqrt(2*units.eV*Ux_in_ExB_frame/units.mp)*1e-3;
  vy_in_ExB = sqrt(2*units.eV*Uy_in_ExB_frame/units.mp)*1e-3;
  vz_in_ExB = sqrt(2*units.eV*Uz_in_ExB_frame/units.mp)*1e-3;
  
  [X,Y,Z] = sphere(200);
  hold(hca,'on')
  %hsurf = surf(hca,vExB_mean(1) + vx_in_ExB*X,vExB_mean(2) + vy_in_ExB*Y, vExB_mean(3) + vz_in_ExB*Z,...
    %'facealpha',0.5,'facecolor',[0 0 0],'EdgeColor','none','FaceLighting','flat');
  hsurf = surf(hca,vExB_mean(1) + vabs_ExB*X,vExB_mean(2) + vabs_ExB*Y, vExB_mean(3) + vabs_ExB*Z,...
    'facealpha',0.5,'facecolor',[0 0.2 0],'EdgeColor','none','FaceLighting','flat');
  %hsurf = surf(hca,vExB_mean(1) + vExB_mean(1)*X,vExB_mean(2) + vExB_mean(2)*Y, vExB_mean(3) + vExB_mean(3)*Z,'facealpha',0.2,'facecolor',[0 0 1],'EdgeColor','none');
  
  hold(hca,'off')

  if 1 % Add uncertainty/std for vExB
    
    [X,Y,Z] = sphere(100);
    hold(hca,'on')
    hsurf = surf(hca,vExB_mean(1) + vExB_std(1)*X,vExB_mean(2) + vExB_std(2)*Y, vExB_mean(3) + vExB_std(3)*Z,...
      'facealpha',0.4,'facecolor',[1 0 0],'EdgeColor','none');
    hold(hca,'off')
  end
  hca.Box = 'on';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  %view(hca,[0 0 1])
  view(hca,[1 1 0.5])
end

if 0
  % create artificial 3D spgerical function that has the value 'iso_values'
  % at a distance defined by the constant surface energy in the ExB frame.
  fff = @(x,y,z,x0,y0,z0,scale,fref) (1/exp(-1))*fref*exp(-(x-x0).^2/scale^2 - (y-y0).^2/scale^2 - (z-z0).^2/scale^2);

  x = linspace(-15000,15000,200);
  y = linspace(-15000,15000,201);
  z = linspace(-15000,15000,202);
  [X,Y,Z] = meshgrid(x,y,z);

  X = hs.grid.vx*1;
  Y = hs.grid.vy*1;
  Z = hs.grid.vz*1;

  f_sphere = fff(X,Y,Z,vExB_mean(1),vExB_mean(2),vExB_mean(3),sqrt(sum(vExB_mean.^2)),iso_values);
  

  nX = permute(X,[2 1 3]);
  nY = permute(Y,[2 1 3]);
  nZ = permute(Z,[2 1 3]);
  % Find the difference field.
  df = hs.data - f_sphere;

  [ii] = find(df < prctile(df,5));
  hold(hca,'on')
  plot3(hca,X(ii),Y(ii),Z(ii),'k.')
  hold(hca,'off')
  %%
  % Interpolate the difference field on the explicitly defined surface.
  f3s = interp3(X, Y, Z, df, X, Y, Z);
  % Find the contour where the difference (on the surface) is zero.
  C = contours(X, Y, Z, f3s, [0 0 0]);
  % Extract the x- and y-locations from the contour matrix C.
  xL = C(1, 2:end);
  yL = C(2, 2:end);
  % Interpolate on the first surface to find z-locations for the intersection
  % line.
  zL = interp2(x2, y2, z2, xL, yL);
  % Visualize the line.
  line(xL,yL,zL,'Color','k','LineWidth',3);
  
  hold(hca,'on')
  isosurf = isosurface(X, Y, Z, f_sphere, iso_values);
  patch(hca,isosurf, 'FaceColor', [0.5 0.0 0.5], 'EdgeColor', 'none');
  hold(hca,'off')
  %%
  patch(isosurface(x3, y3, z3, f2, 0), 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');

  % Define the input grid (in 3D)
  [x3, y3, z3] = meshgrid(linspace(-1,1));
  % Compute the implicitly defined function x^2 + y^2 + z^3 - 0.5^2 = 0
  f1 = x3.^2 + y3.^2 + z3.^2 - 0.5^2;
  % This surface is z = 2*y - 6*x^3, which can also be expressed as
  % 2*y - 6*x^3 - z = 0.
  f2 = 2*y3 - 6*x3.^3 - z3;
  % Also compute z = 2*y - 6*x^3 in the 'traditional' way.
  [x2, y2] = meshgrid(linspace(-1,1));
  z2 = 2*y2 - 6*x2.^3;
  % Visualize the two surfaces.
  patch(isosurface(x3, y3, z3, f1, 0), 'FaceColor', [0.5 1.0 0.5], 'EdgeColor', 'none');
  patch(isosurface(x3, y3, z3, f2, 0), 'FaceColor', [1.0 0.5 0.0], 'EdgeColor', 'none');
  view(3); camlight; axis vis3d;
  % Find the difference field.
  f3 = f1 - f2;
  % Interpolate the difference field on the explicitly defined surface.
  f3s = interp3(x3, y3, z3, f3, x2, y2, z2);
  % Find the contour where the difference (on the surface) is zero.
  C = contours(x2, y2, f3s, [0 0]);
  % Extract the x- and y-locations from the contour matrix C.
  xL = C(1, 2:end);
  yL = C(2, 2:end);
  % Interpolate on the first surface to find z-locations for the intersection
  % line.
  zL = interp2(x2, y2, z2, xL, yL);
  % Visualize the line.
  line(xL,yL,zL,'Color','k','LineWidth',3);
end
if 0 % Add max energy
  %%
  Emax = max(pdist_all.depend{1}(:));
  %Emax = 32000;
  vmax = sqrt(2*units.eV*Emax/units.mp)*1e-3;
  
  [X,Y,Z] = sphere(20);
  hold(hca,'on')
  hsurf = surf(hca,vmax*X,vmax*Y,vmax*Z,'facealpha',0.1,'facecolor',[0 0 0],'EdgeColor','none');
  hold(hca,'off')
end


%cn.print(sprintf([str_tint '_3D_pickup_dt=%03.0fs'],dt))

%% Illustrate motion on a sphere

hca = subplot(1,1,1);
[X,Y,Z] = sphere(200);
v0 = [0.5 -1 0];
v0abs = sqrt(sum(v0.^2));
hsurf = surf(hca,v0(1) + v0abs*X, v0(2) + v0abs*Y, v0(3) + v0abs*Z,...
      'facealpha',0.5,'facecolor',[0.8 0 0.2],'EdgeColor','none','FaceLighting','flat');

axis(hca,'equal')
axis(hca,'square')
camlight(hca,0,0)
   
hca.XLim = 2*[-1 1];
hca.YLim = 2*[-1 1];
hca.ZLim = 2*[-1 1];

hca.XLabel.String = 'v_x';
hca.YLabel.String = 'v_y';
hca.ZLabel.String = 'v_z';

hca.Box = 'on';