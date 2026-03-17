steinThresh = 0.4;
% 1) distances (temperature/covariance similarity)
sd = stein_distance(gm(:,iK), 'sort', true);
% 2) overlap groups (transitive)
groups = overlap_groups_from_distance(sd.D, steinThresh);
% 3) merge groups (moment-matched)
merged = gmm_merge_components(gm(:,iK), groups, 'isort', sd.isort);
% 4) identify cold on merged GMM
coldMerged = gmm_identify_cold_components(merged.gmMerged, ...
  'method','knee','gapPick','last','minGapLog10',0.3,'maxColdK',2,'minColdK',0, ...
  'minWeight',0.01,'minVabs',300);
% 5) keep only cold components that were NOT merged (singleton groups)
cold = gmm_cold_singletons_from_merged(merged, coldMerged);
% cold.coldOriginal{it} and cold.isColdOriginal_mat are what you want

isCold = irf.ts_scalar(times,cellfun(@(x) double(~isempty(x)),coldMerged.coldMerged_kept));

%tsCold_knee = irf.ts_scalar(times,int64(outCold_knee.isCold_mat_sorted));
%tsCold_knee.data(tsCold_knee.data==0) = NaN;

%tsCold_thresh = irf.ts_scalar(times,int64(outCold_thresh.isCold_mat_sorted));
%tsCold_thresh.data(tsCold_thresh.data==0) = NaN;
% Cold indices in sorted component space (what you use with 'sort',true elsewhere)
%cold_sorted_at_t1 = outCold.cold_sorted{1};

% Choose which cold identification to use for masking the overplotted lines
tsCold_mask = tsCold;
%tsCold_mask = tsCold_thresh;

tsVcold = cell(1,K);
for ik = 1:K
  tsVcold{ik} = tsV{ik};
  isNotCold = find(tsCold_mask.data(:,ik)==0);
  tsVcold{ik}.data(isNotCold,:)= NaN;
end