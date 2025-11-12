ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z');

c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)


%% Remove noise and plot
PD_orig = iPDist3;
nMean = [3,3,3,3]; nThresh = 2;
nMean = [5,3,3,3]; nThresh = 3;
nMean = [7,3,3,3]; nThresh = 4;
PD_clean = PD_orig.remove_noise([4 3 3 3],4,iPDist3_counts);
%PD_clean = iPDist3.remove_noise([7 1 1 1],2,iPDist3_counts);
PD_clean = PD_orig.remove_noise(nMean,nThresh,iPDist3_counts);
PD_diff = PD_orig+-PD_clean;

tsElow = PD_orig.find_noise_energy_limit(5).movmean(30);

emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit
n_orig = PD_orig.mask({emask_mat}).n;
n_clean = PD_clean.mask({emask_mat}).n;
n_diff = PD_diff.mask({emask_mat}).n;

v_orig = PD_orig.mask({emask_mat}).vel;
v_clean = PD_clean.mask({emask_mat}).vel;
v_diff = PD_diff.mask({emask_mat}).vel;

%[tsEigVal_orig,tsEigV1_orig,tsEigV2_orig] = PD_orig.movmean(7).p.eig([1 2]); % based on xy pressure
%[tsEigVal_clean,tsEigV1_clean,tsEigV2_clean] = PD_clean.movmean(7).p.eig([1 2]); % based on xy pressure


%
h = irf_plot(5);
fontsize = 14;

if 1 % deflux ion original
  hca = irf_panel('ion orig omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_orig.deflux.nan2zero.omni.specrec,'donotfitcolorbarlabel');
  irf_legend(hca,'Original',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  hca.YScale = 'log';
end
if 1 % deflux ion
  hca = irf_panel('ion counts clean omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_new.deflux.omni.specrec,'donotfitcolorbarlabel');
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Cleaned',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  hca.YScale = 'log';
end
if 1 % deflux ion
  hca = irf_panel('ion counts diff omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_diff.deflux.omni.specrec,'donotfitcolorbarlabel');
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Difference',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  hca.YScale = 'log';
end

if 1 % compare densities
  hca = irf_panel('densities');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{ne3,n_orig,n_clean,n_diff},'comp')
  hca.YLabel.String = 'n (cc)';
  irf_legend(hca,{'ne (FPI)','   Original','   Cleaned','   Difference'},[0.06 0.98],'fontsize',fontsize-2)
end
if 1 % compare vx, less lines
  hca = irf_panel('vx');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{dbcsVi3.x,v_orig.x,v_clean.x,v_diff.x},'comp')
  hca.YLabel.String = 'v_x (km/s)';
  irf_legend(hca,{'v_i (FPI)','   Original','   Cleaned','   Difference'},[0.06 0.98],'fontsize',fontsize-2)
end
if 0 % compare pressures
  hca = irf_panel('x');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_plot(hca,{dbcsVi3.x,v_orig.x,v_clean.x,v_diff.x},'comp')
  hca.YLabel.String = 'v_x (km/s)';
  irf_legend(hca,{'v_i (FPI)','   Original','   Cleaned','   Difference'},[0.06 0.98],'fontsize',fontsize-2)
end

legs = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)'};
for ip = 1:numel(h); irf_legend(h(ip),legs{ip},[-0.08 0.98],'color','k','fontsize',fontsize); end
%h = irf_plot({PD_orig.deflux.omni.specrec,PD_clean.deflux.omni.specrec,PD_diff.deflux.omni.specrec});
irf_plot_axis_align
irf_zoom(h,'x',[PD_orig.time.start PD_orig.time.stop])
h(end).XTickLabelRotation = 0;
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).Layer = ''top'';',1:numel(h))
colormap(irf_colormap('magma'))
hlinks = linkprop(h(1:3),{'CLim','YLim','XLim'});
