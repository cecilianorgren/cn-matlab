% Recalculate moments, all the MP's are saved within moms
iK = 1;
mm = gmm_moments_from_macroparticles(allMP,gm(:,iK),{1,2,3,4,5},'sort',true);
ik = 1;
[tsN,tsV,tsT,tsFx,tsFy,tsFz,tsN_tot,tsFx_tot,tsFy_tot,tsFz_tot] = gmm_get_moments_and_dists(times,gm(:,iK),Fobs.mid{:},ntot,'sort',true);

mpN = irf.ts_scalar(times,squeeze(cellfun(@(x) x.n,mm(:,1,ik))));

%mm_ = mm;%(:,1,ik);
%mpN = irf.ts_scalar(times,cellfun(@(x) x.n,mm_));
h = irf_plot(1);
hca = irf_panel('n');
hca.ColorOrder = mms_colors('1234ab');
irf_plot({tsN{ik},mpN},'comp')

%%
iK = 3;
[tsN,tsV,tsT,tsFx,tsFy,tsFz,tsN_tot,tsFx_tot,tsFy_tot,tsFz_tot] = gmm_get_moments_and_dists(times,gm(:,iK),Fobs.mid{:},ntot,'sort',true);

K = vecK(iK);
minGapLog10 = 0.3;
outCold_knee = gmm_identify_cold_components(gm(:,iK), 'method','knee', 'minGapLog10',minGapLog10, 'minWeight',0.01, 'maxColdK',2,'minColdK',0,'minVabs',300);
outCold_thresh = gmm_identify_cold_components(gm(:,iK), 'method','threshold', 'maxcoldev',1000, 'minWeight',0.01, 'maxColdK',2,'minColdK',0,'minVabs',300);

tsCold_knee = irf.ts_scalar(times,int64(outCold_knee.isCold_mat_sorted));
%tsCold_knee.data(tsCold_knee.data==0) = NaN;

tsCold_thresh = irf.ts_scalar(times,int64(outCold_thresh.isCold_mat_sorted));
%tsCold_thresh.data(tsCold_thresh.data==0) = NaN;
% Cold indices in sorted component space (what you use with 'sort',true elsewhere)
%cold_sorted_at_t1 = outCold.cold_sorted{1};

% Choose which cold identification to use for masking the overplotted lines
tsCold_mask = tsCold_knee;
%tsCold_mask = tsCold_thresh;

tsVcold = cell(1,K);
for ik = 1:K
  tsVcold{ik} = tsV{ik};
  isNotCold = find(tsCold_mask.data(:,ik)==0);
  tsVcold{ik}.data(isNotCold,:)= NaN;
end

ts_gaps = irf.ts_scalar(times,cat(1,[outCold_knee.gaps_log10{:}])');

% Cold indices in the *original* unsorted gm component numbering
%cold_unsorted_at_t1 = outCold.cold_unsorted{1};
%%
figure(13)
h = irf_plot(4);


hca = irf_panel('fred Fz tot');
specrec = tsFz_tot.specrec('velocity');
irf_spectrogram(hca,specrec,'log')    
hca.YLabel.String = 'v_z (km/s)';
irf_legend(hca,'GMM',[0.98 0.98],'k')
hca.NextPlot = "add";
ts = cellfun(@(x){x.z},tsVcold);
hh = irf_plot(hca,ts,'comp');
hca.NextPlot = "replace";

hca = irf_panel('T sorted');
irf_plot(hca,cellfun(@(x){x.trace/3},tsT),'comp');

hca = irf_panel('log gaps');
irf_plot(hca,ts_gaps);
hca.NextPlot = "add";
hh = irf_plot(hca,irf.ts_scalar(times,zeros(times.length,1)+minGapLog10),'k--');
hca.NextPlot = "replace";


hca = irf_panel('v abs');
irf_plot(hca,cellfun(@(x){x.abs},tsV),'comp');

if 0
hca = irf_panel('is cold knee');
irf_plot(hca,tsCold_knee,'*');
end
if 0
hca = irf_panel('is cold thresh');
irf_plot(hca,tsCold_thresh,'*');
end

irf_plot_axis_align(h)
%%

moms = cell(nt,nK);
k = [3:K];
for it = 1:nt
  %MP = moms{it,iK};
  % This can be done afterwards, because R is the same for all K.
  idx = gm{it,iK}.cluster(R);
  MP.("cluster"){K} = idx;
  moms_{it,iK} = gmm_moments_from_macroparticles(MP,idx,1);
end

mpN = irf.ts_scalar(times,cellfun(@(x) x.n,moms));
mpJ = irf.ts_vec_xyz(times,[cellfun(@(x) x.jx,moms) cellfun(@(x) x.jy,moms) cellfun(@(x) x.jz,moms)]);
%mpP = irf.ts_tensor_xyz(times,[cellfun(@(x) jx.n,moms) cellfun(@(x) x.jy,moms) cellfun(@(x) x.jz,moms)]);

h = irf_plot(5);

hca = irf_plot(n);

%% Check the grouping
bd = bhattacharyya_distance(gm, 'sort', true); % sort by gmm_sort (for now scalar temperature)
groups = overlap_groups_from_distance(bd.DB, 0.7);
merged = gmm_merge_components(gm(:,iK), groups, 'isort', bd.isort);
gm_merged = merged.gmMerged;   % cell array (nt x 1) of merged gmdistribution

it = 50;


%%
h = irf_plot(4);
hca = irf_panel('T sorted');
irf_plot(hca,cellfun(@(x){x.trace/3},tsT),'comp');

%%