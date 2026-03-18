steinThresh = 0.3;
steinThresh = 0.3;
% 1) distances (temperature/covariance similarity)
sd = stein_distance(gm, 'sort', true);
% 2) overlap groups (transitive)
groups = overlap_groups_from_distance(sd.D, steinThresh);
% 3) merge groups (moment-matched)
merged = gmm_merge_components(gm, groups, 'isort', sd.isort);
% 4) identify cold on merged GMM
coldMerged = gmm_identify_cold_components(merged.gmMerged, ...
  'method','knee','gapPick','last','minGapLog10',0.3,'maxColdK',2,'minColdK',0, ...
  'minWeight',0.01,'minVabs',300);
% 5) keep only cold components that were NOT merged (singleton groups)
cold = gmm_cold_singletons_from_merged(merged, coldMerged);
% cold.coldOriginal{it} and cold.isColdOriginal_mat are what you want

isCold = irf.ts_scalar(times,cellfun(@(x) double(~isempty(x)),cold.coldMerged_kept));

%tsCold_knee = irf.ts_scalar(times,int64(outCold_knee.isCold_mat_sorted));
%tsCold_knee.data(tsCold_knee.data==0) = NaN;

%tsCold_thresh = irf.ts_scalar(times,int64(outCold_thresh.isCold_mat_sorted));
%tsCold_thresh.data(tsCold_thresh.data==0) = NaN;
% Cold indices in sorted component space (what you use with 'sort',true elsewhere)
%cold_sorted_at_t1 = outCold.cold_sorted{1};

% Choose which cold identification to use for masking the overplotted lines
%tsCold_mask = tsCold;
%tsCold_mask = tsCold_thresh;

%tsVcold = cell(1,K);
%for ik = 1:K
%  tsVcold{ik} = tsV{ik};
%  isNotCold = find(tsCold_mask.data(:,ik)==0);
%  tsVcold{ik}.data(isNotCold,:)= NaN;
%end

%% Make a time loop to plot the results
directory_ = strrep(printpath,'\','');

[h1,h2] = initialize_combined_plot('topbottom',2,1,4,0.5,'vertical');

if 1 % B
  hca = irf_panel('B');
  hca.ColorOrder = mms_colors('xyza');
  irf_plot(hca,{gseB.x, gseB.y, gseB.z},'comp')
  hca.YLabel.String = 'B (nT)';
  hca.ColorOrder = mms_colors('xyza');
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
end
if 1 % fred by mms
  hca = irf_panel('fred z mms');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  irf_spectrogram(hca,vdf_fz.specrec,'log')   
  hca.YLabel.String = 'v_z (km/s)';  
  irf_legend(hca,'MMS',[0.98 0.98],'k')
end

irf_zoom(h1,'x',tint)
irf_plot_axis_align(h1)
h1(end).XTickLabelRotation = 0;

% Plot distributions and GMM results

iK = 3;
K = vecK(iK);

for it = 1:5:nt
  if exist('hmark','var'); delete(hmark); end
  c_eval('hmark = irf_pl_mark(h1,times(it),[0.5 0.5 0.5]);',1:numel(h1))

  isub = 1;
  gm_orig_tmp = gm{it,iK};
  gm_merged_tmp = merged.gmMerged{it,iK};

  hca = h2(isub); isub = isub + 1;
  plot(hca,vdf_fz.depend{1}(it,:),vdf_fz.data(it,:))
  hca.Title.String = 'Observed';
  hca.XLabel.String = 'v (km/s)';
  hca.YLabel.String = sprintf('f (%s)',vdf_fz.units);

  hca = h2(isub); isub = isub + 1;  
  dv = 50;
  vvec = -2500:50:2500;
  [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
  gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
  gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
  plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
  hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
  hca.XLabel.String = 'v (km/s)';
  hca.YLabel.String = sprintf('f (%s)','...');
  legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
  irf_legend(hca,legs',[0.98 0.98])

  hca = h2(isub); isub = isub + 1;  
  dv = 50;
  vvec = -2500:50:2500;
  [gmmFtot, gmmFcomp] = gmm_get_F(gm_merged_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
  gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
  gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
  plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
  hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
  hca.XLabel.String = 'v (km/s)';
  hca.YLabel.String = sprintf('f (%s)','...');
  legs = cellfun(@(x) sprintf('%g',x),groups{it,iK},'UniformOutput',false);
  irf_legend(gca,legs',[0.98 0.98])


  hca = h2(isub); isub = isub + 1;  
  dv = 50;
  vvec = -2500:50:2500;
  [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
  nGroups = numel(groups{it,iK});
  gmmFcomp_z = zeros(nGroups,numel(vvec));
  for iGroup = 1:nGroups
    
  end
  gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
  gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
  plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
  hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
  hca.XLabel.String = 'v (km/s)';
  hca.YLabel.String = sprintf('f (%s)','...');
  legs = cellfun(@(x) sprintf('%g',x),groups{it,iK},'UniformOutput',false);
  irf_legend(gca,legs',[0.98 0.98])

  c_eval('axis(h2(?),''square'');',1:numel(h2))
  cn.print(sprintf('gmm_iDF=%04.f_it=%04.f_K=%g_merged_stein_thresh=%.2f',iDF,it,K,steinThresh))
end