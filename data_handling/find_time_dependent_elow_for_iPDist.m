%% 

%% Plot omni distributions of DEFlux, DPflux, counts

%h = setup_subplots(3,1);

omni_psd = iPDist3.omni;
omni_deflux = iPDist3.deflux.omni;
omni_counts = iPDist3_counts.omni*4*pi; omni_counts.name = 'counts'; omni_counts,units = '';
omni_counts2 = iPDist3_counts.omni; omni_counts2.data = nansum(iPDist3_counts.data(:,:,:),3); omni_counts2.name = 'counts'; omni_counts2,units = '';

nMovMean = 5;
omni_counts_movmean = omni_counts; 
omni_counts_movmean.data = movmean(omni_counts_movmean.data,nMovMean,1);

omni_counts2.data(omni_counts2.data==0) = NaN;
omni_counts2_movmean = omni_counts2; 
omni_counts2_movmean.data = movmean(omni_counts2_movmean.data,nMovMean,1);


%irf_plot({omni_psd.specrec,omni_deflux.specrec,omni_counts.specrec})


if 0 % based on .omni
  log_nan = isnan(omni_counts_movmean.data);
  idx_nan = zeros(omni_counts_movmean.length,1);
  elow_t = zeros(omni_counts_movmean.length,1);
  for it = 1:omni_counts_movmean.length
    idx_nan(it) = find(log_nan(it,:),1,'last');
    elow_t(it) = omni_counts_movmean.depend{1}(it,idx_nan(it));
  end
  nSmoothElow = 10;
  tsElow = irf.ts_scalar(omni_counts_movmean.time,smooth(elow_t,nSmoothElow)); tsElow.name = 'E_low';
end
if 1 % based on sum of counts
  log_nan = isnan(omni_counts_movmean.data);
  idx_nan = zeros(omni_counts_movmean.length,1);
  elow_t = zeros(omni_counts_movmean.length,1);
  for it = 1:omni_counts_movmean.length
    idx_nan(it) = find(log_nan(it,:),1,'last');
    elow_t(it) = omni_counts_movmean.depend{1}(it,idx_nan(it));
  end
  nSmoothElow = 30;
  tsElow = irf.ts_scalar(omni_counts_movmean.time,smooth(elow_t,nSmoothElow)); tsElow.name = 'E_low';
end
if 0 % based on sum of counts
  min_counts = 10;
  data = omni_counts_movmean.data;
  elow_t = zeros(omni_counts_movmean.length,1);
  clear idx_
  for it = 1:omni_counts_movmean.length
    idx_(it) = find(data(it,:) > min_counts,1,'first');
    if isempty(idx_(it)), idx_(it) = size(data,2); end

    elow_t(it) = omni_counts_movmean.depend{1}(it,idx_(it));
  end
  nSmoothElow = 10;
  tsElow = irf.ts_scalar(omni_counts_movmean.time,smooth(elow_t,nSmoothElow)); tsElow.name = 'E_low';
end


%last_idx_nan = find(idx_nan);

fontsize = 12;
h = irf_plot(4);

if 1
  hca = irf_panel('omni psd');
  irf_spectrogram(hca,omni_psd.specrec,'log','donotfitcolorbarlabel')
  hca.YScale = 'log';
end
if 1
  hca = irf_panel('omni deflux');
  irf_spectrogram(hca,omni_deflux.specrec,'log','donotfitcolorbarlabel')
  hca.YScale = 'log';
end
if 0
  hca = irf_panel('omni counts');
  [~,hcb] = irf_spectrogram(hca,omni_counts.specrec,'lin','donotfitcolorbarlabel');
  hca.YScale = 'log';
  hca.CLim(1) = 0;
end
if 1
  hca = irf_panel('omni counts 2');
  [~,hcb] = irf_spectrogram(hca,omni_counts2.specrec,'lin','donotfitcolorbarlabel');
  hca.YScale = 'log';
  hca.CLim(1) = 0;
  hcb.YLabel.String = 'Total counts';
end

if 1
  hca = irf_panel('omni counts movmean 5');
  [~,hcb] = irf_spectrogram(hca,omni_counts2_movmean.specrec,'lin','donotfitcolorbarlabel');
  hca.YScale = 'log';
  hcb.YLabel.String = 'Total counts';
  
  irf_legend(hca,{sprintf('%g-point averaged counts',nMovMean),sprintf('%g-point averaged E_{low}',nSmoothElow)}',[0.02 0.1],'color','k','fontsize',fontsize,'fontweight','bold','backgroundcolor','w')

  hca.CLim(1) = 0;
  hold(hca,'on')
  irf_plot(tsElow,'k','linewidth',2)
  hold(hca,'off')
  hca.YLabel.String = 'E_{low} (eV)';
  hca.YLabel.Interpreter = 'tex';
end

c_eval('h(?).Layer = "top";',1:numel(h))
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
h(end).XTickLabelRotation = 0;
h(1).Title.String = 'Finding time-variable lower energy limit';