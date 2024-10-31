%%
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT') +- 0;
nMovMean = 7;
elim = [3000 Inf];
pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).elim(elim).tlim(time_xline + 0.15*0.5*[-1 1]);
ipdist = iPDist3.movmean(nMovMean).tlim(time_xline + 0.15*0.5*[-1 1]);
%epdist = ePDist3.movmean(nMovMean).tlim(time_xline + 0.03*0.5*[-1 1]);

%%
N = 100000;
MP = pdist.elim([50 Inf]).macroparticles('ntot',N,'skipzero',1);
V = [MP.vx, MP.vy, MP.vz];

%
nGroups = 4;  % number of classes
[idx, c, sumd, d] = kmeans(V, nGroups);
%s = silhouette(V, idx, 'sqeuclid');

%%
idx = clusterR;
fieldnames = fields(MP);
for iGroup = 1:nGroups
  idx_tmp = find(idx==iGroup); 
  for iField = 1:numel(fieldnames)
    var = MP.(fieldnames{iField});
    MPtmp.(fieldnames{iField}) = var(idx_tmp);
  end
  
  % Calculate speed and density
  n = sum(MPtmp.dn);
  jx = sum(MPtmp.vx.*MPtmp.dn);
  jy = sum(MPtmp.vy.*MPtmp.dn);
  jz = sum(MPtmp.vz.*MPtmp.dn);
  vx = jx/n;
  vy = jy/n;
  vz = jz/n;
  MPtmp.sum_n = n*1e-6;
  MPtmp.sum_jx = jx;
  MPtmp.sum_jy = jy;
  MPtmp.sum_jz = jz;
  MPtmp.sum_vx = vx;
  MPtmp.sum_vy = vy;
  MPtmp.sum_vz = vz;

  MP_grouped{iGroup} = MPtmp;
end

%% Make nGroups pdists based on the dominant class for each bin.
Dep1_edges = [pdist.ancillary.energy-pdist.ancillary.delta_energy_minus pdist.ancillary.energy(end)+pdist.ancillary.delta_energy_minus(end)];
Dep2_edges = 0.5:1:32.5;
Dep3_edges = 0.5:1:16.5;
for iGroup = 1:nGroups
  [average_groupid,~,mid,~] = histcn([MP.iDep1, MP.iDep2, MP.iDep3], Dep1_edges, Dep2_edges, Dep3_edges, 'AccumData', idx, 'Fun', @mean);
  average_groupid = round(average_groupid); 
  pdist_tmp = pdist;
  id_exclude = find(average_groupid ~= iGroup);
  %[i1,i2,i3] = ind2sub([32 32 16],id_exclude);
  pdist_tmp.data(1,id_exclude) = 0;
  
  pdist_group{iGroup} = pdist_tmp;
end


% klist = 2:4; 
% myfunc = @(X,K)(kmeans(X, K));
% eva = evalclusters(net.IW{1},myfunc,'CalinskiHarabasz','klist',klist);
% classes = kmeans(net.IW{1},eva.OptimalK);
% 
% eva = evalclusters(V,'kmeans','CalinskiHarabasz','KList',1:6);
% Optimal_K = eva.OptimalK;

%% 
h = setup_subplots(2,2);
isub = 1;

hca = h(isub); isub = isub + 1;
hs = scatter3(hca, V(:,1), V(:,2), V(:,3), 15, idx);
hca.XLabel.String = 'v_x (km/s)';
hca.YLabel.String = 'v_y (km/s)';
hca.ZLabel.String = 'v_z (km/s)';
cmap = colormap;
icolors = round(interp1(1:size(cmap,1),1:size(cmap,1),linspace(1,size(cmap,1),nGroups)));
group_colors = cmap(icolors,:);


if 1 % f(x,y)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',[1 0 0],[0 1 0]);
  vdf.plot_plane(hca);
  hca.XLabel.String = 'v_x (km/s)';
  hca.YLabel.String = 'v_y (km/s)';
  hold(hca,'on')
  for iGroup = 1:nGroups
    plot(hca,MP_grouped{iGroup}.sum_vx,MP_grouped{iGroup}.sum_vy,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
  end
  hold(hca,'off')
end
if 1 % f(x,z)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',[1 0 0],[0 0 1]);
  vdf.plot_plane(hca);
  hca.XLabel.String = 'v_x (km/s)';
  hca.YLabel.String = 'v_z (km/s)';
  hold(hca,'on')
  for iGroup = 1:nGroups
    plot(hca,MP_grouped{iGroup}.sum_vx,MP_grouped{iGroup}.sum_vz,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
  end
  hold(hca,'off')
end
if 1 % f(y,z)
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',[0 1 0],[0 0 1]);
  vdf.plot_plane(hca);
  hca.XLabel.String = 'v_y (km/s)';
  hca.YLabel.String = 'v_z (km/s)';
  hold(hca,'on')
  for iGroup = 1:nGroups
    plot(hca,MP_grouped{iGroup}.sum_vy,MP_grouped{iGroup}.sum_vz,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
  end
  hold(hca,'off')
end

%%
dv = 100; % cm/s
nscaling = 1/(0.01^5);
for iGroup = 1:nGroups
  hca = h(isub); isub = isub + 1;
  v_edges = -2200:dv:2200;
  MPtmp = MP_grouped{iGroup};
  [count,~,mid,~] = histcn([MPtmp.vx, MPtmp.vy], v_edges, v_edges, 'AccumData', MPtmp.df*nscaling, 'Fun', @sum);
  count(count==0) = NaN;
  pcolor(hca,mid{1},mid{2},log10(count*dv*1e2)')
  hca.XLabel.String = 'v_x (km/s)';
  hca.YLabel.String = 'v_y (km/s)';
  shading(hca,'flat')
  hcb = colorbar(hca);
end

for iGroup = 1:nGroups
  hca = h(isub); isub = isub + 1;
  v_edges = -2200:dv:2200;
  MPtmp = MP_grouped{iGroup};
  [count,~,mid,~] = histcn([MPtmp.vx, MPtmp.vz], v_edges, v_edges, 'AccumData', MPtmp.df*nscaling, 'Fun', @sum);
  count(count==0) = NaN;
  pcolor(hca,mid{1},mid{2},log10(count*dv*1e2)')
  hca.XLabel.String = 'v_x (km/s)';
  hca.YLabel.String = 'v_z (km/s)';
  shading(hca,'flat')
  hcb = colorbar(hca);
end

for iGroup = 1:nGroups
  hca = h(isub); isub = isub + 1;
  v_edges = -2200:dv:2200;
  MPtmp = MP_grouped{iGroup};
  [count,~,mid,~] = histcn([MPtmp.vy, MPtmp.vz], v_edges, v_edges, 'AccumData', MPtmp.df*nscaling, 'Fun', @sum);
  count(count==0) = NaN;
  pcolor(hca,mid{1},mid{2},log10(count*dv*1e2)')
  hca.XLabel.String = 'v_y (km/s)';
  hca.YLabel.String = 'v_z (km/s)';
  shading(hca,'flat')
  hcb = colorbar(hca);
end


%% Run through several times
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT') +- 0;
tint_kmeans = time_xline + 10*[-1 1];

nMovMean = 7;
elim = [3000 Inf];
pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).elim(elim).tlim(tint_kmeans);

nGroups = 4;
nMacroParticles = 20000;


initial_guess = [0 -1000 -1000;...
                 0 -1000 +1000;...
                 0 +1000 -1000;...
                 0 +1000 +1000];

times = EpochTT([]);
tic
for iPD = 1:pdist.length
  pdist_tmp = pdist(iPD);
  kmp = kmeans_pdist(pdist_tmp,nGroups,nMacroParticles,'start', initial_guess);
  
  times = [times kmp.t];
  for iGroup = 1:nGroups
    n(iPD,iGroup) = kmp.MP_grouped{iGroup}.sum_n;
    vx(iPD,iGroup) = kmp.MP_grouped{iGroup}.sum_vx;
    vy(iPD,iGroup) = kmp.MP_grouped{iGroup}.sum_vy;
    vz(iPD,iGroup) = kmp.MP_grouped{iGroup}.sum_vz;
  end  
end
toc

%%
tsN = irf.ts_scalar(pdist.time,n);
tsVx = irf.ts_scalar(pdist.time,vx);
tsVy = irf.ts_scalar(pdist.time,vy);
tsVz = irf.ts_scalar(pdist.time,vz);

tsJx = irf.ts_scalar(pdist.time,n.*vx);
tsJy = irf.ts_scalar(pdist.time,n.*vy);
tsJz = irf.ts_scalar(pdist.time,n.*vz);

n_sum_mat = repmat(sum(n,2),[1 nGroups]);
tsVx_rel = irf.ts_scalar(pdist.time,n.*vx./n_sum_mat);
tsVy_rel = irf.ts_scalar(pdist.time,n.*vy./n_sum_mat);
tsVz_rel = irf.ts_scalar(pdist.time,n.*vz./n_sum_mat);

%% Plot results
h = irf_plot(7);

hca = irf_panel('n');
irf_plot(hca,{tsN},'*');
hca.YLabel.String = 'n (cc?)';

if 1
  hca = irf_panel('vx');
  irf_plot(hca,{tsVx},'*');
  hca.YLabel.String = 'v_x (km/s)';
end
if 1
  hca = irf_panel('vy');
  irf_plot(hca,{tsVy},'*');
  hca.YLabel.String = 'v_y (km/s)';
end
if 1
  hca = irf_panel('vz');
  irf_plot(hca,{tsVz},'*');
  hca.YLabel.String = 'v_z (km/s)';
end
if 1
  hca = irf_panel('jx');
  irf_plot(hca,{(tsJx)},'*');
  hca.YLabel.String = 'j_x (cc*km/s)';
end
if 1
  hca = irf_panel('jy');
  irf_plot(hca,{tsJy},'*');
  hca.YLabel.String = 'j_y (cc*km/s)';
end
if 1
  hca = irf_panel('jz');
  irf_plot(hca,{tsJz},'*');
  hca.YLabel.String = 'j_z (cc*km/s)';
end

irf_pl_mark(h,time_xline,'k')
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
