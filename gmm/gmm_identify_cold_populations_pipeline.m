
flag_distance = 'stein';               distanceThresh = 0.2;
flag_distance = 'bhattacharyya';       distanceThresh = 0.4;
flag_distance = 'bhattacharyya_term1'; distanceThresh = 0.3;
flag_distance = 'gaussianity';         distanceThresh = 0.3;

% 1) distances (temperature/covariance similarity)
switch flag_distance
  case 'stein'
    sd = stein_distance(gm, 'sort', true);
  case 'bhattacharyya'
    sd = bhattacharyya_distance(gm, 'sort', true);
  case 'bhattacharyya_term1'
    sd = onlydrift_distance(gm, 'sort', true);
  case 'gaussianity'
    tic;
    sd = merge_gaussians_maxwellianity(gm,'method','grid'); 
    t_mc = toc;
end
%%

% ..) find how the maxwellianity changes with different levels of merging
%     two ways to do this, iteratively, with iterative merging, or checking
%     all piece-wise measures first, and then find merging points from
%     that (from no merges to all merging)

successive_groups = successive_groups_from_distance(sd.D);
successive_merged = gmm_merge_components(gm, successive_groups); % merge by law of shared covariances etc...
tic;
successive_quality = gmm_compare_vdfs_envelope(gm,successive_merged.gmMerged,'method_points','grid');
toc
threshold_quality = 0.1;
% find threshold for which the "quality" is below a given value



distanceThresh = 0.2;
% 2a) find the number of final goups based on threshold
distanceThreshVec = 0:0.1:1;
ngroups_vs_threshold = ngroups_from_threshold(sd.D,distanceThreshVec);
% 2) overlap groups (transitive)
groups = overlap_groups_from_distance(sd.D, distanceThresh);
% 3) merge groups (moment-matched)
%merged = gmm_merge_components(gm, groups, 'isort', sd.isort);
merged = gmm_merge_components(gm, groups); % merge by law of shared covariances etc...
% 4) identify cold on merged GMM
coldMerged = gmm_identify_cold_components(merged.gmMerged, ...
  'method','knee','gapPick','last','minGapLog10',10^-0.5,'maxColdK',2,'minColdK',0, ...
  'minWeight',0.01,'minVabs',200);
% 5) keep only cold components that were NOT merged (singleton groups)
cold = gmm_cold_singletons_from_merged(merged, coldMerged);
% cold.coldOriginal{it} and cold.isColdOriginal_mat are what you want

isCold = irf.ts_scalar(times,cellfun(@(x) double(~isempty(x)),cold.coldMerged_kept));


% Calculate the final difference between the original gm and the merged gm



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

%%  Make a time loop to plot the results
directory_ = strrep(printpath,'\','');

iK = 1;
K = vecK(iK);

ts_ngroups_vs_threshold = irf.ts_scalar(times,cat(1,ngroups_vs_threshold{:,iK}));
ts_ngroups_vs_threshold.userData.specrec_f = distanceThreshVec;

ts_grouping_quality = irf.ts_scalar(times,squeeze(successive_quality(:,iK,:)));


[h1,h2] = initialize_combined_plot('topbottom',5,3,4,0.5,'horizontal');
fontsize = 10;

if 1 % B
  hca = irf_panel('B');
  hca.ColorOrder = mms_colors('xyza');
  irf_plot(hca,{gseB.x, gseB.y, gseB.z},'comp')
  hca.YLabel.String = 'B (nT)';
  hca.ColorOrder = mms_colors('xyza');
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98],'fontsize',fontsize-1)
end
if 1 % fred by mms
  hca = irf_panel('fred z mms');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  irf_spectrogram(hca,vdf_fz.specrec,'log')   
  hca.YLabel.String = 'v_z (km/s)';  
  irf_legend(hca,'MMS',[0.98 0.98],'color','k','fontsize',fontsize-1)
end
if 1 % fred by mms
  hca = irf_panel('is cold');
  hca.ColorOrder = pic_colors('matlab');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  irf_plot(hca,isCold,'*')
  hca.YLabel.String = 'cold flag';  
  legs = "K=" + cellfun(@(x)x.NumComponents,gm(1,:));
  irf_legend(hca,legs,[0.98 0.98],'color','k','fontsize',fontsize-1)
end
if 1 % n groups based on threshold
  hca = irf_panel('ngroups based on threshold');
  hca.ColorOrder = pic_colors('matlab');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  [hax,hcb] = irf_spectrogram(hca,ts_ngroups_vs_threshold.specrec,'lin');
  hcb.YLim = [1 K];
  %legs = "K=" + cellfun(@(x)x.NumComponents,gm(1,:));
  %irf_legend(hca,legs,[0.98 0.98],'k')
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,ones(times.length,1)*distanceThresh),'k')
  hold(hca,'off')
  hca.YLabel.String = {'Distance','threshold',flag_distance,sprintf('K=%g',K)};
end
if 1 % n groups based on threshold
  hca = irf_panel('quality based on n groups');
  hca.ColorOrder = pic_colors('matlab');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  [hax,hcb] = irf_spectrogram(hca,ts_grouping_quality.specrec,'lin');
  hca.YLim = [1 K]+[-0.5 0.5];
  hcb.Label.String = 'RMS diff';
  %legs = "K=" + cellfun(@(x)x.NumComponents,gm(1,:));
  %irf_legend(hca,legs,[0.98 0.98],'k')
  %hold(hca,'on')
  %irf_plot(hca,irf.ts_scalar(times,ones(times.length,1)*distanceThresh),'k')
  %hold(hca,'off')
  hca.YLabel.String = {'N_{groups}',sprintf('K=%g',K)}';
end


%irf_zoom(h1,'x',tint)
irf_plot_axis_align(h1)
h1(end).XTickLabelRotation = 0;

% Plot distributions and GMM results
iG = 1;
for it = 20;81;%1:5:nt
  if exist('hmark','var'); delete(hmark); end
  c_eval('hmark = irf_pl_mark(h1,times(it),[0.5 0.5 0.5]);',1:numel(h1))

  isub = 1;
  gm_orig_tmp = gm{it,iK};
  gm_merg_tmp = merged.gmMerged{it,iK};


  dv = 50;
  vvec = -2500:50:2500;
  [gmmFtot, gmmFcomp, gmmFgrouped] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it),'group',groups{it,iK}{iG});  
  gmmFtot_x = squeeze(sum(gmmFtot,[2 3]))*(dv*dv*1e6)*1e-18;
  gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;
  gmmFgrouped_x = squeeze(sum(gmmFgrouped,[2 3]))*dv*dv*1e6*1e-18;
  gmmFtot_y = squeeze(sum(gmmFtot,[1 3]))*(dv*dv*1e6)*1e-18;
  gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;
  gmmFgrouped_y = squeeze(sum(gmmFgrouped,[1 3]))*dv*dv*1e6*1e-18;
  gmmFtot_z = squeeze(sum(gmmFtot,[1 2]))*(dv*dv*1e6)*1e-18;
  gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;
  gmmFgrouped_z = squeeze(sum(gmmFgrouped,[1 2]))*dv*dv*1e6*1e-18;

  [gmmFtot_merg, gmmFcomp_merg] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it));
  gmmFcomp_x_merg = squeeze(sum(gmmFcomp_merg,[2 3]))*(dv*dv*1e6)*1e-18;
  gmmFcomp_y_merg = squeeze(sum(gmmFcomp_merg,[1 3]))*(dv*dv*1e6)*1e-18;
  gmmFcomp_z_merg = squeeze(sum(gmmFcomp_merg,[1 2]))*(dv*dv*1e6)*1e-18;

  gmmFtot_x_merg = squeeze(sum(gmmFtot_merg,[2 3]))*(dv*dv*1e6)*1e-18;
  gmmFtot_y_merg = squeeze(sum(gmmFtot_merg,[1 3]))*(dv*dv*1e6)*1e-18;
  gmmFtot_z_merg = squeeze(sum(gmmFtot_merg,[1 2]))*(dv*dv*1e6)*1e-18;

  dd = gmm_compare_vdfs(gm_orig_tmp,gm_merg_tmp);

  if 1 % f(vx)
    % Original gmm
    dv = 50;
    vvec = -2500:50:2500;
    %[gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    %gmmFtot_x = squeeze(sum(gmmFtot,[2 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
    %gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
  
    hca = h2(isub); isub = isub + 1;
    hca.ColorOrder = mms_colors('1234');
    plot(hca,vdf_fx.depend{1}(it,:),vdf_fx.data(it,:),vvec,gmmFtot_x,vvec,gmmFtot_x_merg)
    hca.Title.String = 'Observed + Total GMM';
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = sprintf('f (%s)',vdf_fx.units);
    hca.ColorOrder = mms_colors('1234');
    irf_legend(hca,{'Obs.','GMM','GMM merged',sprintf('diff = %5.3f',dd)}',[0.02 0.98],'fontsize',fontsize-1)
  
    hca = h2(isub); isub = isub + 1;  
    plot(hca,vvec,gmmFcomp_x,vvec,gmmFtot_x,'k')
    hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    if 1 % print the distances
      %toprint = cellstr("" + sd.D{it,iK});
      toprint = arrayfun(@(x) sprintf('%01.2f',x),sd.D{it,iK},'UniformOutput',false);
      toprint = num2str(sd.D{it,iK}, '%5.2f');
      irf_legend(hca,toprint,[0.02 0.98],'color','k','fontsize',9)    
    end
 
    hca = h2(isub); isub = isub + 1;  
    % dv = 50;
    % vvec = -2500:dv:2500;
    % [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    % gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;
    % nGroups = numel(groups{it,iK});
    % gmmFcomp_x_grouped = zeros(numel(vvec),nGroups);
    % for iGroup = 1:nGroups
    %   group_tmp = groups{it,iK}{iGroup};
    %   gmmFcomp_x_grouped(:,iGroup) = sum(gmmFcomp_x(:,group_tmp),2);    
    % end
    plot(hca,vvec,gmmFgrouped_x,vvec,sum(gmmFgrouped_x,2),'k',vdf_fx.depend{1}(it,:),vdf_fx.data(it,:),'k--')
    %gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
    %gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
    %plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
    %
    hca.Title.String = {'Summed f',sprintf('D_{%s}<%g',flag_distance,distanceThresh)};
    %hca.Title.String = {sprintf('Gaussian Mixture Model'),'summed components'};
    hca.XLabel.String = 'v (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    %irf_legend(hca,{'Summed','components'}',[0.02 0.98],'color','k')

    hca = h2(isub); isub = isub + 1;  
    dv = 50;
    vvec = -2500:50:2500;
    %[gmmFtot, gmmFcomp] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    %gmmFtot_x = squeeze(sum(gmmFtot,[2 3]))*(dv*dv*1e6)*1e-18;
    %gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;
    plot(hca,vvec,gmmFcomp_x_merg,vvec,gmmFtot_x_merg,'k',vdf_fx.depend{1}(it,:),vdf_fx.data(it,:),'k--')
    %hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
    hca.Title.String = {'Merged by law of shared \mu, \Sigma, w',sprintf('D_{Stein}<%g',distanceThresh)};
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    %irf_legend(hca,{'Merged by','law of','merged','covariances'}',[0.02 0.98],'color','k')
  

  end

  if 1 % f(vy)
    % Original gmm
    % dv = 50;
    % vvec = -2500:50:2500;
    % [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    % gmmFtot_y = squeeze(sum(gmmFtot,[1 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
    % gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
    % 
    hca = h2(isub); isub = isub + 1;
    hca.ColorOrder = mms_colors('1234');
    plot(hca,vdf_fy.depend{1}(it,:),vdf_fy.data(it,:),vvec,gmmFtot_y,vvec,gmmFtot_y_merg)
    hca.Title.String = 'Observed + Total GMM';
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = sprintf('f (%s)',vdf_fy.units);
    hca.ColorOrder = mms_colors('1234');
    irf_legend(hca,{'Obs.','GMM','GMM merged'}',[0.02 0.98],'fontsize',fontsize-1)
  
    hca = h2(isub); isub = isub + 1;  
    plot(hca,vvec,gmmFcomp_y,vvec,gmmFtot_y,'k')
    hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    if 1 % print the distances
      %toprint = cellstr("" + sd.D{it,iK});
      toprint = arrayfun(@(x) sprintf('%01.2f',x),sd.D{it,iK},'UniformOutput',false);
      toprint = num2str(sd.D{it,iK}, '%5.2f');
      irf_legend(hca,toprint,[0.02 0.98],'color','k','fontsize',9)    
    end

    hca = h2(isub); isub = isub + 1;  
    % dv = 50;
    % vvec = -2500:dv:2500;
    % [gmmFtot, gmmFcomp, gmmFgrouped] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it),'group',groups{it,iK}); 
    % gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;
    % nGroups = numel(groups{it,iK});
    % gmmFcomp_y_grouped = zeros(numel(vvec),nGroups);
    % for iGroup = 1:nGroups
    %   group_tmp = groups{it,iK}{iGroup};
    %   gmmFcomp_y_grouped(:,iGroup) = sum(gmmFcomp_y(:,group_tmp),2);    
    % end
    % gmmFcomp_y_grouped = squeeze(sum(gmmFgrouped,[1 3]))*dv*dv*1e6*1e-18;
    plot(hca,vvec,gmmFgrouped_y,vvec,sum(gmmFgrouped_y,2),'k',vdf_fy.depend{1}(it,:),vdf_fy.data(it,:),'k--')
    %gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
    %gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
    %plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
    %
    hca.Title.String = {'Summed f',sprintf('D_{%s}<%g',flag_distance,distanceThresh)};
    %hca.Title.String = {sprintf('Gaussian Mixture Model'),'summed components'};
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    %irf_legend(hca,{'Summed','components'}',[0.02 0.98],'color','k')

    hca = h2(isub); isub = isub + 1;  
    % dv = 50;
    % vvec = -2500:50:2500;
    % [gmmFtot, gmmFcomp] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    % gmmFtot_y = squeeze(sum(gmmFtot,[1 3]))*(dv*dv*1e6)*1e-18;
    % gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;
    plot(hca,vvec,gmmFcomp_y_merg,vvec,gmmFtot_y_merg,'k',vdf_fy.depend{1}(it,:),vdf_fy.data(it,:),'k--')
    %hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
    hca.Title.String = {'Merged by law of shared \mu, \Sigma, w',sprintf('D_{Stein}<%g',distanceThresh)};
    hca.XLabel.String = 'v (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    %irf_legend(hca,{'Merged by','law of','merged','covariances'}',[0.02 0.98],'color','k')
  

  end


  if 1 % f(vz)
    % Original gmm
    dv = 50;
    vvec = -2500:50:2500;
    [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    gmmFtot_z = squeeze(sum(gmmFtot,[1 2]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
    gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
  
    hca = h2(isub); isub = isub + 1;
    hca.ColorOrder = mms_colors('1234');
    plot(hca,vdf_fz.depend{1}(it,:),vdf_fz.data(it,:),vvec,gmmFtot_z,vvec,gmmFtot_z_merg)
    hca.Title.String = 'Observed + Total GMM';
    hca.XLabel.String = 'v_z (km/s)';
    hca.YLabel.String = sprintf('f (%s)',vdf_fz.units);
    hca.ColorOrder = mms_colors('1234');
    irf_legend(hca,{'Obs.','GMM','GMM merged'}',[0.02 0.98],'fontsize',fontsize-1)
  
    hca = h2(isub); isub = isub + 1;  
    plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
    hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
    hca.XLabel.String = 'v_z (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    if 1 % print the distances
      %toprint = cellstr("" + sd.D{it,iK});
      toprint = arrayfun(@(x) sprintf('%01.2f',x),sd.D{it,iK},'UniformOutput',false);
      toprint = num2str(sd.D{it,iK}, '%5.2f');
      irf_legend(hca,toprint,[0.02 0.98],'color','k','fontsize',9)    
    end
  
    hca = h2(isub); isub = isub + 1;  
    % dv = 50;
    % vvec = -2500:dv:2500;
    % [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    % gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;
    % nGroups = numel(groups{it,iK});
    % gmmFcomp_z_grouped = zeros(numel(vvec),nGroups);
    % for iGroup = 1:nGroups
    %   group_tmp = groups{it,iK}{iGroup};
    %   gmmFcomp_z_grouped(:,iGroup) = sum(gmmFcomp_z(:,group_tmp),2);    
    % end
    plot(hca,vvec,gmmFgrouped_z,vvec,sum(gmmFgrouped_z,2),'k',vdf_fz.depend{1}(it,:),vdf_fz.data(it,:),'k--')
    %gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
    %gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
    %plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
    %
    hca.Title.String = {'Summed f',sprintf('D_{%s}<%g',flag_distance,distanceThresh)};
    %hca.Title.String = {sprintf('Gaussian Mixture Model'),'summed components'};
    hca.XLabel.String = 'v (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    %irf_legend(hca,{'Summed','components'}',[0.02 0.98],'color','k')

    hca = h2(isub); isub = isub + 1;  
    % dv = 50;
    % vvec = -2500:50:2500;
    % [gmmFtot, gmmFcomp] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
    % gmmFtot_z = squeeze(sum(gmmFtot,[1 2]))*(dv*dv*1e6)*1e-18;
    % gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;
    plot(hca,vvec,gmmFcomp_z_merg,vvec,gmmFtot_z_merg,'k',vdf_fz.depend{1}(it,:),vdf_fz.data(it,:),'k--')
    %hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
    hca.Title.String = {'Merged by law of shared \mu, \Sigma, w',sprintf('D_{Stein}<%g',distanceThresh)};
    hca.XLabel.String = 'v_z (km/s)';
    hca.YLabel.String = sprintf('f (%s)','...');
    legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
    irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
    %irf_legend(hca,{'Merged by','law of','merged','covariances'}',[0.02 0.98],'color','k')
  

  end

  c_eval('h1(?).FontSize = fontsize;',1:numel(h1));
  c_eval('h2(?).FontSize = fontsize;',1:numel(h2));
  c_eval('h2(?).Title = [];',5:numel(h2));
  c_eval('axis(h2(?),''square'');',1:numel(h2))
  %cn.print(sprintf('gmm_iDF=%04.f_it=%04.f_K=%g_merged_stein_thresh=%.2f',iDF,it,K,steinThresh))
end


%% Other type of loop, to illustrate successive merging
% Make a time loop to plot the results
directory_ = strrep(printpath,'\','');

iK = 1;
K = vecK(iK);

ts_ngroups_vs_threshold = irf.ts_scalar(times,cat(1,ngroups_vs_threshold{:,iK}));
ts_ngroups_vs_threshold.userData.specrec_f = distanceThreshVec;

ts_grouping_quality = irf.ts_scalar(times,squeeze(successive_quality(:,iK,:)));


%[h1,h2] = initialize_combined_plot('topbottom',3,5,4,0.3,'horizontal');
[h1,h2] = initialize_combined_plot('leftright',3,5,3,0.4,'horizontal');
fontsize = 10;

if 1 % B
  hca = irf_panel('B');
  hca.ColorOrder = mms_colors('xyza');
  irf_plot(hca,{gseB.x, gseB.y, gseB.z},'comp')
  hca.YLabel.String = 'B (nT)';
  hca.ColorOrder = mms_colors('xyza');
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98],'fontsize',fontsize-1)
end
if 1 % fred by mms
  hca = irf_panel('fred z mms');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  irf_spectrogram(hca,vdf_fz.specrec,'log')   
  hca.YLabel.String = 'v_z (km/s)';  
  irf_legend(hca,'MMS',[0.98 0.98],'color','k','fontsize',fontsize-1)
end
if 0 % fred by mms
  hca = irf_panel('is cold');
  hca.ColorOrder = pic_colors('matlab');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  irf_plot(hca,isCold,'*')
  hca.YLabel.String = 'cold flag';  
  legs = "K=" + cellfun(@(x)x.NumComponents,gm(1,:));
  irf_legend(hca,legs,[0.98 0.98],'color','k','fontsize',fontsize-1)
end
if 0 % n groups based on threshold
  hca = irf_panel('ngroups based on threshold');
  hca.ColorOrder = pic_colors('matlab');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  [hax,hcb] = irf_spectrogram(hca,ts_ngroups_vs_threshold.specrec,'lin');
  hcb.YLim = [1 K];
  %legs = "K=" + cellfun(@(x)x.NumComponents,gm(1,:));
  %irf_legend(hca,legs,[0.98 0.98],'k')
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,ones(times.length,1)*distanceThresh),'k')
  hold(hca,'off')
  hca.YLabel.String = {'Distance','threshold',flag_distance,sprintf('K=%g',K)};
end
if 1 % n groups based on threshold
  hca = irf_panel('quality based on n groups');
  hca.ColorOrder = pic_colors('matlab');
  %vdfx = PD.reduce('1D',[1 0 0]);
  %vdfx = PD.reduce('1D',[1 0 0]);
  % reduce just one before the loops
  [hax,hcb] = irf_spectrogram(hca,ts_grouping_quality.specrec,'lin');
  hca.YLim = [1 K];
  hcb.Label.String = 'RMS diff';
  %legs = "K=" + cellfun(@(x)x.NumComponents,gm(1,:));
  %irf_legend(hca,legs,[0.98 0.98],'k')
  %hold(hca,'on')
  %irf_plot(hca,irf.ts_scalar(times,ones(times.length,1)*distanceThresh),'k')
  %hold(hca,'off')
  hca.YLabel.String = {'N_{groups}',sprintf('K=%g',K)}';
end


irf_zoom(h1,'x',tint)
irf_plot_axis_align(h1)
h1(end).XTickLabelRotation = 0;

% Plot distributions and GMM results
iG = 1;
for it = 20;81;%1:5:nt
  if exist('hmark','var'); delete(hmark); end
  c_eval('hmark = irf_pl_mark(h1,times(it),[0.5 0.5 0.5]);',1:numel(h1))

  isub = 1;
  gm_orig_tmp = gm{it,iK};

  for iG = 1:K
    gm_merg_tmp = successive_merged.gmMerged{it,iK,iG};
  
  
    dv = 50;
    vvec = -2500:50:2500;
    [gmmFtot, gmmFcomp, gmmFgrouped] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it),'group',successive_groups{it,iK}{iG});  
    gmmFtot_x = squeeze(sum(gmmFtot,[2 3]))*(dv*dv*1e6)*1e-18;
    gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;
    gmmFgrouped_x = squeeze(sum(gmmFgrouped,[2 3]))*dv*dv*1e6*1e-18;
    gmmFtot_y = squeeze(sum(gmmFtot,[1 3]))*(dv*dv*1e6)*1e-18;
    gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;
    gmmFgrouped_y = squeeze(sum(gmmFgrouped,[1 3]))*dv*dv*1e6*1e-18;
    gmmFtot_z = squeeze(sum(gmmFtot,[1 2]))*(dv*dv*1e6)*1e-18;
    gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;
    gmmFgrouped_z = squeeze(sum(gmmFgrouped,[1 2]))*dv*dv*1e6*1e-18;
  
    [gmmFtot_merg, gmmFcomp_merg] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it));
    gmmFcomp_x_merg = squeeze(sum(gmmFcomp_merg,[2 3]))*(dv*dv*1e6)*1e-18;
    gmmFcomp_y_merg = squeeze(sum(gmmFcomp_merg,[1 3]))*(dv*dv*1e6)*1e-18;
    gmmFcomp_z_merg = squeeze(sum(gmmFcomp_merg,[1 2]))*(dv*dv*1e6)*1e-18;
  
    gmmFtot_x_merg = squeeze(sum(gmmFtot_merg,[2 3]))*(dv*dv*1e6)*1e-18;
    gmmFtot_y_merg = squeeze(sum(gmmFtot_merg,[1 3]))*(dv*dv*1e6)*1e-18;
    gmmFtot_z_merg = squeeze(sum(gmmFtot_merg,[1 2]))*(dv*dv*1e6)*1e-18;
  
    dd = gmm_compare_vdfs(gm_orig_tmp,gm_merg_tmp);
  
    if 0 % f(vx)
      % Original gmm
      dv = 50;
      vvec = -2500:50:2500;
      %[gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      %gmmFtot_x = squeeze(sum(gmmFtot,[2 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
      %gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
    
      hca = h2(isub); isub = isub + 1;
      hca.ColorOrder = mms_colors('1234');
      plot(hca,vdf_fx.depend{1}(it,:),vdf_fx.data(it,:),vvec,gmmFtot_x,vvec,gmmFtot_x_merg)
      hca.Title.String = 'Observed + Total GMM';
      hca.XLabel.String = 'v_x (km/s)';
      hca.YLabel.String = sprintf('f (%s)',vdf_fx.units);
      hca.ColorOrder = mms_colors('1234');
      irf_legend(hca,{'Obs.','GMM','GMM merged',sprintf('diff = %5.3f',dd)}',[0.02 0.98],'fontsize',fontsize-1)
    
      hca = h2(isub); isub = isub + 1;  
      plot(hca,vvec,gmmFcomp_x,vvec,gmmFtot_x,'k')
      hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
      hca.XLabel.String = 'v_x (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      if 1 % print the distances
        %toprint = cellstr("" + sd.D{it,iK});
        toprint = arrayfun(@(x) sprintf('%01.2f',x),sd.D{it,iK},'UniformOutput',false);
        toprint = num2str(sd.D{it,iK}, '%5.2f');
        irf_legend(hca,toprint,[0.02 0.98],'color','k','fontsize',9)    
      end
   
      hca = h2(isub); isub = isub + 1;  
      % dv = 50;
      % vvec = -2500:dv:2500;
      % [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      % gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;
      % nGroups = numel(groups{it,iK});
      % gmmFcomp_x_grouped = zeros(numel(vvec),nGroups);
      % for iGroup = 1:nGroups
      %   group_tmp = groups{it,iK}{iGroup};
      %   gmmFcomp_x_grouped(:,iGroup) = sum(gmmFcomp_x(:,group_tmp),2);    
      % end
      plot(hca,vvec,gmmFgrouped_x,vvec,sum(gmmFgrouped_x,2),'k',vdf_fx.depend{1}(it,:),vdf_fx.data(it,:),'k--')
      %gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
      %gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
      %plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
      %
      hca.Title.String = {'Summed f',sprintf('D_{%s}<%g',flag_distance,distanceThresh)};
      %hca.Title.String = {sprintf('Gaussian Mixture Model'),'summed components'};
      hca.XLabel.String = 'v (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      %irf_legend(hca,{'Summed','components'}',[0.02 0.98],'color','k')
  
      hca = h2(isub); isub = isub + 1;  
      dv = 50;
      vvec = -2500:50:2500;
      %[gmmFtot, gmmFcomp] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      %gmmFtot_x = squeeze(sum(gmmFtot,[2 3]))*(dv*dv*1e6)*1e-18;
      %gmmFcomp_x = squeeze(sum(gmmFcomp,[2 3]))*(dv*dv*1e6)*1e-18;
      plot(hca,vvec,gmmFcomp_x_merg,vvec,gmmFtot_x_merg,'k',vdf_fx.depend{1}(it,:),vdf_fx.data(it,:),'k--')
      %hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
      hca.Title.String = {'Merged by law of shared \mu, \Sigma, w',sprintf('D_{Stein}<%g',distanceThresh)};
      hca.XLabel.String = 'v_x (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      %irf_legend(hca,{'Merged by','law of','merged','covariances'}',[0.02 0.98],'color','k')
    
  
    end
  
    if 0 % f(vy)
      % Original gmm
      % dv = 50;
      % vvec = -2500:50:2500;
      % [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      % gmmFtot_y = squeeze(sum(gmmFtot,[1 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
      % gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
      % 
      hca = h2(isub); isub = isub + 1;
      hca.ColorOrder = mms_colors('1234');
      plot(hca,vdf_fy.depend{1}(it,:),vdf_fy.data(it,:),vvec,gmmFtot_y,vvec,gmmFtot_y_merg)
      hca.Title.String = 'Observed + Total GMM';
      hca.XLabel.String = 'v_y (km/s)';
      hca.YLabel.String = sprintf('f (%s)',vdf_fy.units);
      hca.ColorOrder = mms_colors('1234');
      irf_legend(hca,{'Obs.','GMM','GMM merged'}',[0.02 0.98],'fontsize',fontsize-1)
    
      hca = h2(isub); isub = isub + 1;  
      plot(hca,vvec,gmmFcomp_y,vvec,gmmFtot_y,'k')
      hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
      hca.XLabel.String = 'v_y (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      if 1 % print the distances
        %toprint = cellstr("" + sd.D{it,iK});
        toprint = arrayfun(@(x) sprintf('%01.2f',x),sd.D{it,iK},'UniformOutput',false);
        toprint = num2str(sd.D{it,iK}, '%5.2f');
        irf_legend(hca,toprint,[0.02 0.98],'color','k','fontsize',9)    
      end
  
      hca = h2(isub); isub = isub + 1;  
      % dv = 50;
      % vvec = -2500:dv:2500;
      % [gmmFtot, gmmFcomp, gmmFgrouped] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it),'group',groups{it,iK}); 
      % gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;
      % nGroups = numel(groups{it,iK});
      % gmmFcomp_y_grouped = zeros(numel(vvec),nGroups);
      % for iGroup = 1:nGroups
      %   group_tmp = groups{it,iK}{iGroup};
      %   gmmFcomp_y_grouped(:,iGroup) = sum(gmmFcomp_y(:,group_tmp),2);    
      % end
      % gmmFcomp_y_grouped = squeeze(sum(gmmFgrouped,[1 3]))*dv*dv*1e6*1e-18;
      plot(hca,vvec,gmmFgrouped_y,vvec,sum(gmmFgrouped_y,2),'k',vdf_fy.depend{1}(it,:),vdf_fy.data(it,:),'k--')
      %gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
      %gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
      %plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
      %
      hca.Title.String = {'Summed f',sprintf('D_{%s}<%g',flag_distance,distanceThresh)};
      %hca.Title.String = {sprintf('Gaussian Mixture Model'),'summed components'};
      hca.XLabel.String = 'v_y (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      %irf_legend(hca,{'Summed','components'}',[0.02 0.98],'color','k')
  
      hca = h2(isub); isub = isub + 1;  
      % dv = 50;
      % vvec = -2500:50:2500;
      % [gmmFtot, gmmFcomp] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      % gmmFtot_y = squeeze(sum(gmmFtot,[1 3]))*(dv*dv*1e6)*1e-18;
      % gmmFcomp_y = squeeze(sum(gmmFcomp,[1 3]))*(dv*dv*1e6)*1e-18;
      plot(hca,vvec,gmmFcomp_y_merg,vvec,gmmFtot_y_merg,'k',vdf_fy.depend{1}(it,:),vdf_fy.data(it,:),'k--')
      %hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
      hca.Title.String = {'Merged by law of shared \mu, \Sigma, w',sprintf('D_{Stein}<%g',distanceThresh)};
      hca.XLabel.String = 'v (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = cellfun(@(x) sprintf('%g',x),groups{it,iK}{iG},'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      %irf_legend(hca,{'Merged by','law of','merged','covariances'}',[0.02 0.98],'color','k')
    
  
    end
  
    if 1 % f(vz)
      % Original gmm
      dv = 50;
      vvec = -2500:50:2500;
      [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      gmmFtot_z = squeeze(sum(gmmFtot,[1 2]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
      gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;%*(dv*dv*1e6)*1e-18;
    
      hca = h2(isub); isub = isub + 1;
      hca.ColorOrder = mms_colors('1234');
      plot(hca,vdf_fz.depend{1}(it,:),vdf_fz.data(it,:),vvec,gmmFtot_z,vvec,gmmFtot_z_merg)
      hca.Title.String = 'Observed + Total GMM';
      hca.XLabel.String = 'v_z (km/s)';
      hca.YLabel.String = sprintf('f (%s)',vdf_fz.units);
      hca.ColorOrder = mms_colors('1234');
      irf_legend(hca,{'Obs.','GMM','GMM merged'}',[0.02 0.98],'fontsize',fontsize-1)
    
      if 0 % original with distnce measure
        hca = h2(isub); isub = isub + 1;  
        plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
        hca.Title.String = {sprintf('Gaussian Mixture Model, K=%g',K)};
        hca.XLabel.String = 'v_z (km/s)';
        hca.YLabel.String = sprintf('f (%s)','...');
        legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
        irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
        if 1 % print the distances
          %toprint = cellstr("" + sd.D{it,iK});
          toprint = arrayfun(@(x) sprintf('%01.2f',x),sd.D{it,iK},'UniformOutput',false);
          toprint = num2str(sd.D{it,iK}, '%5.2f');
          irf_legend(hca,toprint,[0.02 0.98],'color','k','fontsize',9)    
        end
      end
    
      hca = h2(isub); isub = isub + 1;  
      % dv = 50;
      % vvec = -2500:dv:2500;
      % [gmmFtot, gmmFcomp] = gmm_get_F(gm_orig_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      % gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;
      % nGroups = numel(groups{it,iK});
      % gmmFcomp_z_grouped = zeros(numel(vvec),nGroups);
      % for iGroup = 1:nGroups
      %   group_tmp = groups{it,iK}{iGroup};
      %   gmmFcomp_z_grouped(:,iGroup) = sum(gmmFcomp_z(:,group_tmp),2);    
      % end
      plot(hca,vvec,gmmFgrouped_z,vvec,sum(gmmFgrouped_z,2),'k',vdf_fz.depend{1}(it,:),vdf_fz.data(it,:),'k--')
      %gmmFtot_z = squeeze(sum(gmmFtot,[1 2]));%*(dv*dv*1e6)*1e-18;
      %gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]));%*(dv*dv*1e6)*1e-18;
      %plot(hca,vvec,gmmFcomp_z,vvec,gmmFtot_z,'k')
      %
      hca.Title.String = {'Summed f'};
      %hca.Title.String = {sprintf('Gaussian Mixture Model'),'summed components'};
      hca.XLabel.String = 'v_z (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = cellfun(@(x) sprintf('%g',x),successive_groups{it,iK}{iG},'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      %irf_legend(hca,{'Summed','components'}',[0.02 0.98],'color','k')
  
      hca = h2(isub); isub = isub + 1;  
      % dv = 50;
      % vvec = -2500:50:2500;
      % [gmmFtot, gmmFcomp] = gmm_get_F(gm_merg_tmp,vvec,vvec,vvec,ntot(it)); %[X,Y,Z] = ndgrid(Fobs.mid{:});
      % gmmFtot_z = squeeze(sum(gmmFtot,[1 2]))*(dv*dv*1e6)*1e-18;
      % gmmFcomp_z = squeeze(sum(gmmFcomp,[1 2]))*(dv*dv*1e6)*1e-18;
      plot(hca,vvec,gmmFcomp_z_merg,vvec,gmmFtot_z_merg,'k',vdf_fz.depend{1}(it,:),vdf_fz.data(it,:),'k--')
      %hca.Title.String = {sprintf('Gaussian Mixture Model'),'merged \mu, \Sigma, w'};
      hca.Title.String = {'Merged by law of shared \mu, \Sigma, w'};
      hca.XLabel.String = 'v_z (km/s)';
      hca.YLabel.String = sprintf('f (%s)','...');
      legs = cellfun(@(x) sprintf('%g',x),successive_groups{it,iK}{iG},'UniformOutput',false);
      irf_legend(hca,legs',[0.98 0.98],'fontsize',fontsize-1)
      %irf_legend(hca,{'Merged by','law of','merged','covariances'}',[0.02 0.98],'color','k')
    
  
    end
  end
  c_eval('h1(?).FontSize = fontsize;',1:numel(h1));
  c_eval('h2(?).FontSize = fontsize;',1:numel(h2));
  c_eval('h2(?).Title = [];',5:numel(h2));
  c_eval('axis(h2(?),''square'');',1:numel(h2))
  %cn.print(sprintf('gmm_iDF=%04.f_it=%04.f_K=%g_merged_stein_thresh=%.2f',iDF,it,K,steinThresh))
end

