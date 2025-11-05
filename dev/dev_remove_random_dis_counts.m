% 
PD = iPDist3;
PD_counts = iPDist3_counts;
elim_lowE = [0 100];
av_counts = mean(sum(PD_counts.nan2zero.data,[3 4]),1);
counts_25 = prctile(sum(PD_counts.nan2zero.data,[3 4]),25,1);
counts_50 = prctile(sum(PD_counts.nan2zero.data,[3 4]),50,1);
counts_75 = prctile(sum(PD_counts.nan2zero.data,[3 4]),75,1);

counts_E = sum(PD_counts.nan2zero.data(:,:,:,:),[3 4]);

av_counts_lowE = mean(sum(PD_counts.elim(elim_lowE).nan2zero.data,[3 4]),1);
av_counts_lowE_1 = mean(av_counts_lowE);

N_rem_per_elevel = round(av_counts_lowE_1);

f_per_count = PD.data./PD_counts.data;

f_per_count_1 = nanmean(nanmean(f_per_count,[3 4]),1);

%%
PD_new = PD;
for it = 1:PD.length
  it
  for iE = 1:32
    if sum(PD_counts.data(it,iE,:)) == N_rem_per_elevel

    nRem = N_rem_per_elevel;
    iRem = 0;
    while iRem < nRem
      iAz = randi(32,1);
      iPol = randi(16,1);
      if PD.data(it,iE,iAz,iPol) == 0
        continue
      else
        PD_new.data(it,iE,iAz,iPol) = PD_new.data(it,iE,iAz,iPol) - f_per_count(it,iE,iAz,iPol);
        iRem = iRem + 1;
      end
    end
  end
end

%%
fontsize = 14;

hca = subplot(2,1,1);
plot(hca,PD_counts.depend{1}(1,:),av_counts)
hold(hca,'on')
hl = plot(hca,PD_counts.elim(elim_lowE).depend{1}(1,:),av_counts_lowE,'*');
hold(hca,'off')
hca.XLabel.String = 'Energy (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = 'Counts summed over angles and averaged over time';
hca.XScale = 'log';
irf_legend(hca,{sprintf('Average at low energies: %.2f',av_counts_lowE_1)},[0.02 0.98],'color',hl.Color,'fontsize',fontsize)
hca.FontSize = fontsize;


hca = subplot(2,1,2);
plot(hca,PD_counts.depend{1}(1,:),[counts_25;counts_50;counts_75])
%hold(hca,'on')
%hl = plot(hca,PD_counts.elim(elim_lowE).depend{1}(1,:),av_counts_lowE,'*');
%hold(hca,'off')
hca.XLabel.String = 'Energy (eV)';
hca.YLabel.String = 'Counts';
hca.Title.String = 'Counts summed over angles and percentiled over time';
hca.XScale = 'log';
%irf_legend(hca,{sprintf('Average at low energies: %.2f',av_counts_lowE_1)},[0.02 0.98],'color',hl.Color,'fontsize',fontsize)
hca.FontSize = fontsize;

%%
fontsize = 14;

hca = subplot(1,1,1);
elim_lowE = [0 100];
counts_E = sum(PD_counts.elim(elim_lowE).nan2zero.data(:,:,:,:),[3 4]);
hl = histogram(hca,counts_E(:),-0.5:1:25,'Normalization','pdf');
irf_legend(hca,{sprintf('%g < E < %g eV',elim_lowE(1),elim_lowE(2))},[0.98 0.98],'fontsize',fontsize)
hca.FontSize = fontsize;
hca.YLabel.String = 'Probability distribution function';
hca.XLabel.String = 'Ion counts';
hca.Title.String = 'Counts per energy level (summed over time)';

%% Smooth data then remove, now in space as well, before it was only averaged in time
nMean = [7 2 2 2];
nMean = [7 2 2 2];
nMean = [5 2 2 2];
%nMean = [4 1 5 5];
nPad = [0 0 nMean(3) 0];
nThresh = 3;

data = PD.data;
counts = PD_counts.nan2zero.data;

data_pad = padarray(data,nPad,'circular');
counts_pad = padarray(counts,nPad,'circular');

data_pad_mean = data_pad;
counts_pad_mean = counts_pad;
counts_pad_sum = counts_pad;
for iDim = 1:numel(nMean)
  data_pad_mean = movmean(data_pad_mean,nMean(iDim),iDim);
  counts_pad_mean = movmean(counts_pad_mean,nMean(iDim),iDim);
  counts_pad_sum = movsum(counts_pad_sum,nMean(iDim),iDim);
  %counts_pad_sum = counts_pad_mean*nMean(iDim);
end
% Remove padded dimensions
data_mean = data_pad_mean(nPad(1)+1:end-nPad(1),nPad(2)+1:end-nPad(2),nPad(3)+1:end-nPad(3),nPad(4)+1:end-nPad(4));
counts_mean = counts_pad_mean(nPad(1)+1:end-nPad(1),nPad(2)+1:end-nPad(2),nPad(3)+1:end-nPad(3),nPad(4)+1:end-nPad(4));
counts_sum = counts_pad_sum(nPad(1)+1:end-nPad(1),nPad(2)+1:end-nPad(2),nPad(3)+1:end-nPad(3),nPad(4)+1:end-nPad(4));

counts_mean_clean = counts_mean;
counts_mean_clean(counts_sum<nThresh) = 0;
counts_clean = counts;
counts_clean(counts_sum<nThresh) = 0;

data_mean_clean = data_mean;
data_mean_clean(counts_sum<nThresh) = 0; 
data_clean = data;
data_clean(counts_sum<nThresh) = 0; 
% Instead of doing this, remove the datapoint that contributed to the
% count from the original PD. Averaging over energies and angles becomes 
% too smoothed. Hm, this will likely become a problem around the adges of
% the real data as well...

PDc2_orig = PD.clone(PD.time,counts);
PDc2_new = PD.clone(PD.time,counts_clean);
PDc2_diff = PD.clone(PD.time,PDc2_orig.data-PDc2_new.data);

PDc_orig = PD.clone(PD.time,counts_mean);
PDc_new = PD.clone(PD.time,counts_mean_clean);
PDc_diff = PD.clone(PD.time,PDc_orig.data-PDc_new.data);

PD_orig = PD.clone(PD.time,data_mean);
PD_new = PD.clone(PD.time,data_mean_clean);
PD_diff = PD.clone(PD.time,PD_orig.data-PD_new.data);

PD2_orig = PD.clone(PD.time,data);
PD2_new = PD.clone(PD.time,data_clean);
PD2_diff = PD.clone(PD.time,PD2_orig.data-PD2_new.data);


fontsize = 18;
h = irf_plot(6);
if 1 % counts ion
  hca = irf_panel('ion counts orig omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_counts.nan2zero.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,'Original',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 1 % counts ion
  hca = irf_panel('ion counts orig averaged omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PDc2_orig.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Averaged',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 1 % counts ion
  hca = irf_panel('ion counts clean omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PDc2_new.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'One-counts removed',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 1 % counts ion
  hca = irf_panel('ion counts diff omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PDc2_diff.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Averaged minus one-counts removed',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 0 % deflux ion
  hca = irf_panel('ion counts orig omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD.deflux.nan2zero.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'DEFlux'};
  irf_legend(hca,'Original',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 0 % deflux ion
  hca = irf_panel('ion counts orig averaged omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_orig.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'DEFlux'};
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Averaged',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 0 % deflux ion
  hca = irf_panel('ion counts clean omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_new.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'DEFlux'};
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'One-counts removed',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 0 % deflux ion
  hca = irf_panel('ion counts diff omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,PD_diff.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'DEFlux'};
  irf_legend(hca,{sprintf('Window = [%.0f,%.0f,%.0f,%.0f]',nMean(1),nMean(2),nMean(3),nMean(4));sprintf('N < %g removed',nThresh)},[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Averaged minus one-counts removed',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end

if 1 % compare densities, averaged
  hca = irf_panel('densities');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{ne3,ni3,PD_orig.elim([1000 Inf]).n,PD_new.elim([1000 Inf]).n,PD_diff.elim([1000 Inf]).n},'comp')
  hca.YLabel.String = 'n (cc)';
  irf_legend(hca,{'ne (FPI)','ni (FPI)','Original','One-counts removed','Averaged minus one-counts removed'},[0.98 0.98])
end
if 1 % compare densities
  hca = irf_panel('densities 2');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{ne3,ni3,PD2_orig.elim([1000 Inf]).n,PD2_new.elim([1000 Inf]).n,PD2_diff.elim([1000 Inf]).n},'comp')
  hca.YLabel.String = 'n (cc)';
  irf_legend(hca,{'ne (FPI)','ni (FPI)','Original','One-counts removed','Averaged minus one-counts removed'},[0.98 0.98])
end
c_eval('h(?).FontSize = fontsize;',1:numel(h))
hlinks = linkprop(h(1:3),{'CLim','YLim'});
colormap(flipdim(irf_colormap('Spectral'),1))
%h(1).CLim = [-4 -0.5];
%h(1).CLim = [-28 -20];
hb = findobj(gcf,'type','colorbar'); %hb = hb(end:-1:1);
%c_eval('hb(?).YLabel.String = ''DEFlux'';',1:numel(hb))
irf_plot_axis_align(h)
irf_zoom(h,'x',[iPDist3.time.start iPDist3.time.stop])
%%

data_pad_mean = smooth3(data_pad,'box',nMean);

data_pad_mean_clean = data_pad_mean;
data_pad_mean_clean(counts_pad_sum<1.1) = NaN;

%% Illustrate what has been removed
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT') +- 0;
nMovMean = 7;
elim = [3000 Inf];

ipdist = iPDist3.movmean(nMovMean);
ipdist_counts = iPDist3_counts.nan2zero.movmean(nMovMean);
ipdist_counts_clean = iPDist3_counts.nan2zero.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts);
ipdist_clean = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts);
ipdist_clean_elim = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).elim(elim);
ipdist_elim = iPDist3.elim(elim);

ipdist_diff = ipdist.clone(ipdist.time,ipdist.data-ipdist_clean.data);
ipdist_counts_diff = ipdist_counts.clone(ipdist_counts.time,ipdist_counts.data-ipdist_counts_clean.data);


%
fontsize = 18;
h = irf_plot(3);

if 1 % counts ion
  hca = irf_panel('ion counts omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,ipdist_counts.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,sprintf('N = %.0f',nMovMean),[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Original',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 1 % counts ion cleaned
  hca = irf_panel('ion counts cleaned omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,ipdist_counts_clean.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,sprintf('N = %.0f',nMovMean),[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'One-counts removed',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end
if 1 % counts ion
  hca = irf_panel('ion counts diff omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  [hs, hcb] = irf_spectrogram(hca,ipdist_counts_diff.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  hcb.YLabel.String = {'Counts'};
  irf_legend(hca,sprintf('N = %.0f',nMovMean),[0.02 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
  irf_legend(hca,'Original minus one-counts removed',[0.98 0.1],'fontsize',fontsize,'color','k','backgroundcolor','w')
end

colormap(flipdim(irf_colormap('Spectral'),1))
c_eval('h(?).FontSize = fontsize;',1:numel(h))
hlinks = linkprop(h,{'CLim','YLim'});
h(1).CLim = [-4 -0.5];