% load /Users/Cecilia/Data/MMS/2015Oct16/workspace2015June15.mat
tint = irf.tint('2015-10-16T12:01:00.00Z/2015-10-16T13:59:00.00Z');
tint = irf.tint('2015-10-16T10:01:00.00Z/2015-10-16T11:59:00.00Z');
ic = 1;

%% Load particle distribution
%c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,ic);',ic)
%c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,ic);',ic)

c_eval('iPDist? = mms.get_data(''PDi_fpi_fast_l2'',tint,ic);',ic)
c_eval('ePDist? = mms.get_data(''PDe_fpi_fast_l2'',tint,ic);',ic)

c_eval('ne? = mms.get_data(''Ne_fpi_fast_l2'',tint,ic);',ic)
c_eval('ni? = mms.get_data(''Ni_fpi_fast_l2'',tint,ic);',ic)

c_eval('vi? = mms.get_data(''Vi_dbcs_fpi_fast_l2'',tint,ic);',ic)
c_eval('ve? = mms.get_data(''Ve_dbcs_fpi_fast_l2'',tint,ic);',ic)

c_eval('Ti? = mms.get_data(''Ti_dbcs_fpi_brst_l2'',tint,ic);',ic)
c_eval('Te? = mms.get_data(''Te_dbcs_fpi_brst_l2'',tint,ic);',ic)

c_eval('Pi? = mms.get_data(''Pi_dbcs_fpi_brst_l2'',tint,ic);',ic)
c_eval('Pe? = mms.get_data(''Pe_dbcs_fpi_brst_l2'',tint,ic);',ic)

load /Users/Cecilia/Data/MMS/2015Oct16/defatt.mat

c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc',ic);
c_eval('tic; gseB?fast=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);

R  = mms.get_data('R_gse',tint);
c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);


%% Make figure
tint = irf.tint('2015-10-16T12:01:00.00Z/2015-10-16T13:59:00.00Z');
tint = irf.tint('2015-10-16T10:01:00.00Z/2015-10-16T11:59:00.00Z');
ic = 1;

h = irf_plot(10);

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?fast.x.tlim(tint),gseB?fast.y.tlim(tint),gseB?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{vi?.x.tlim(tint),vi?.y.tlim(tint),vi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist omni 64
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
%
if 0
if 1 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,18).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
end
if 1 % Ve
  c_eval('Veplot = gseVe?;',ic)
  c_eval('neplot = ne?;',ic)
  tmpdata = Veplot.data;
  indLowNe = find(neplot.data<1);
  tmpdata(indLowNe,:) = repmat([NaN NaN NaN],numel(indLowNe),1);
  Veplot.data =tmpdata;
  c_eval('Veplot = gseVe?.clone(gseVe?.time,tmpdata);',ic)
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{Veplot.x.tlim(tint),Veplot.y.tlim(tint),Veplot.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_e_x','V_e_y','V_e_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'n>1 cm^{-3}'},[0.08 0.95],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [-1200 700];
  
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hout,irf_colormap('space'))  
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist omni 32
  hca = irf_panel('e DEF omni');
  c_eval('irf_spectrogram(hca,ePDist?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.pitchangles(dmpaB?,24).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 0 % i DEF omni 64
  hca = irf_panel('i DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iDEFomni64_?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hout,irf_colormap('space'))  
end
if 0 % i DEF omni 32
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(tint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 0 % J curl  
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
end
%load('caa/cmap.mat');
for ii = [];[4 5 7 8]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  %colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:10
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
end
%colormap('jet')

irf_zoom(h,'x',tint)
%irf_zoom(h([1:3 6 10]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);