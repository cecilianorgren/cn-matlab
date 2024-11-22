%%
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT') +- 0;
nMovMean = 20;
elim = [3000 Inf];
pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).elim(elim).tlim(time_xline + 0.15*0.5*[-1 1]);
ipdist = iPDist3.movmean(nMovMean).tlim(time_xline + 0.15*0.5*[-1 1]);
%epdist = ePDist3.movmean(nMovMean).tlim(time_xline + 0.03*0.5*[-1 1]);

%%
N = 100000;
MP = pdist.elim([50 Inf]).macroparticles('ntot',N,'skipzero',1);
V = [MP.vx, MP.vy, MP.vz];

%
initial_guess = [0 -1000 -1000;...
                 0 -1000 +1000;...
                 0 +1000 -1000;...
                 0 +1000 +1000];

nGroups = 4;  % number of classes
[idx, c, sumd, d] = kmeans(V, nGroups,'Start',initial_guess);
s = silhouette(V, idx, 'sqeuclid');


MP_grouped = macroparticle_moments(MP,idx);
pdist_group = partial_pdist(pdist,MP,idx);

gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],4);

pdist_group_gmm = partial_pdist(pdist,MP,gm_out.idx);
if 0
  %%
  evaluation = evalclusters(V,"kmeans","silhouette","KList",1:8)
end

%% Calculate moments and split macroparticle structure into one for each group
%idx = clusterR;
if 1
  MP_grouped = macroparticle_moments(MP,idx);
else
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
end

%% Make nGroups pdists based on the dominant class for each bin.

if 1
  pdist_group = partial_pdist(pdist,MP,idx);
else
  %Dep1_edges = [pdist.ancillary.energy-pdist.ancillary.delta_energy_minus pdist.ancillary.energy(end)+pdist.ancillary.delta_energy_minus(end)];
  Dep1_edges = 0.5:1:(size(pdist.depend{1},2)+0.5);
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
end


% klist = 2:4; 
% myfunc = @(X,K)(kmeans(X, K));
% eva = evalclusters(net.IW{1},myfunc,'CalinskiHarabasz','klist',klist);
% classes = kmeans(net.IW{1},eva.OptimalK);
% 
% eva = evalclusters(V,'kmeans','CalinskiHarabasz','KList',1:6);
% Optimal_K = eva.OptimalK;

%% Plot
delete(h)
h = setup_subplots(2,1,'vertical');
isub = 1;

colormap(pic_colors('matlab'))

hca = h(isub); isub = isub + 1;
hs = scatter3(hca, V(:,1), V(:,2), V(:,3), 30, idx,'filled');
hca.XLabel.String = 'v_x (km/s)';
hca.YLabel.String = 'v_y (km/s)';
hca.ZLabel.String = 'v_z (km/s)';
hs.MarkerFaceAlpha = 0.1;
hs.MarkerEdgeAlpha = 0.1;

cmap = colormap;
icolors = round(interp1(1:size(cmap,1),1:size(cmap,1),linspace(1,size(cmap,1),nGroups)));
group_colors = cmap(icolors,:);

view(hca,[0 0 1])

%isub = isub + 1;

if 0 % f(x,y)
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
if 0 % f(x,z)
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
if 0 % f(y,z)
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

if 1 % f(x,y), contour plots of all groups together
  hca = h(isub); isub = isub + 1;
  holdon = 0;
  for iGroup = 1:numel(pdist_group)
    %set(hca,'colororder',group_colors(iGroup,:))
    if not(holdon); hold(hca,'on'); holdon = 1; end   
    pd = pdist_group{iGroup};
    vdf = pd.reduce('2D',[1 0 0],[0 1 0]);
    data = squeeze(vdf.data);
    %data(isnan(data)) = 0;
    nSmooth = 0;
    data = smooth2(data,nSmooth,nSmooth);
    %vdf.data(vdf.data==0) = NaN;
    %vdf.data(vdf.data) = NaN;
    data = log10(data);


    [hh_,hh] = contour(vdf.depend{1},vdf.depend{1},data',1,'color',group_colors(iGroup,:),'linewidth',4);
    %[hh_,hh] = contour(vdf.depend{1},vdf.depend{1},data',3,'k');
    %hh.FaceColor = 'none';
    %hh.Color = group_colors(iGroup,:);
    hh.FaceAlpha = 0.2;
    %hh.EdgeColor = group_colors(iGroup,:);
    pause(1)
  end
  hold(hca,'off');
  hca.XLabel.String = 'v_x (km/s)';
  hca.YLabel.String = 'v_y (km/s)';
  hold(hca,'on')
  for iGroup = 1:nGroups
    plot(hca,MP_grouped{iGroup}.sum_vx,MP_grouped{iGroup}.sum_vy,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
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


%% Run through times, do gaussian mixture model and kmeans, and compare


for dt = 7
  time = time_xline + dt;
  nMovMean = 7;
  pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).tlim(time + 0.15*0.5*[-1 1]).elim([300 Inf]);


  N = 100000;
  MP = pdist.elim([50 Inf]).macroparticles('ntot',N,'skipzero',1);
  V = [MP.vx, MP.vy, MP.vz];

  % Initial guess
  initial_guess = [0 -1000 -1000;...
                   0 -1000 +1000;...
                   0 +1000 -1000;...
                   0 +1000 +1000];
  % kmeans  
  nGroups = 4;  % number of classes
  [idx, c, sumd, d] = kmeans(V, nGroups,'Start',initial_guess);
  %s = silhouette(V, idx, 'sqeuclid');
  
  
  MP_grouped_km = macroparticle_moments(MP,idx);
  pdist_group_km = partial_pdist(pdist,MP,idx);
  

  % Gaussian mixture model
  clear S
  S.mu = initial_guess;    
  % S.Sigma = 500*ones(3,3,nGroups);
  S.Sigma = 500*ones(1,3,nGroups); % with input: 'CovType','diagonal'
  S.ComponentProportion = [1 1 1 1];
  
  %gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],4,'Start',S);
  gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],4,'Start',S,'CovType','diagonal');
  
    
  %gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],4);
  
  MP_grouped_gmm = macroparticle_moments(MP,gm_out.idx);
  pdist_group_gmm = partial_pdist(pdist,MP,gm_out.idx);
  

  %% Plot results
  
  nRows = 3;
  nCols = 2;  
  [h1,h] = initialize_combined_plot('topbottom',1,nRows,nCols,0.2,'horizontal');
  h1.Position(2) = h1.Position(2)-0.05;
  h1.Position(4) = h1.Position(4)*2;
  fig = gcf;
  fig.Position = [98         181         751         1078];
  
  hca = h1(1);
  colors = mms_colors('xyza');
  set(hca,'ColorOrder',colors)  
  hh = irf_plot(hca,gseVi3,'linewidth',1);
  c_eval('hh(?).Color = colors(?,:);',1:3)
  hca.YLabel.String = 'v_i (km/s)';
  hca.YLabel.Interpreter = 'tex';
  
  irf_zoom(h1,'x',[gseVi3.time.start gseVi3.time.stop] + 130*[1 -1])
  if exist('hmark','var'); delete(hmark); end
  hmark = irf_pl_mark(h1,pdist.time + nMovMean*0.150*0.5*[-1 1],'k');


  
  colormap(pic_colors('matlab'))
  linewidth_contour = 2;
  vlim = 2500;
  
  isub = 1;
  step = 100; % plot only every xth macroparticles, for speed
  if 1 % Scatter plot of gmm 
    hca = h(isub); isub = isub + 1;
    hs = scatter3(hca, V(1:step:end,1), V(1:step:end,2), V(1:step:end,3), 30, gm_out.idx(1:step:end),'filled');
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    hca.ZLabel.String = 'v_z (km/s)';
    hs.MarkerFaceAlpha = 0.1;
    hs.MarkerEdgeAlpha = 0.1;
    
    cmap = colormap;
    icolors = round(interp1(1:size(cmap,1),1:size(cmap,1),linspace(1,size(cmap,1),nGroups)));
    group_colors = cmap(icolors,:);
    hca.Title.String = 'Gaussian mixture model';
  end
  if 1 % Scatter plot of kmeans 
    hca = h(isub); isub = isub + 1;
    hs = scatter3(hca, V(1:step:end,1), V(1:step:end,2), V(1:step:end,3), 30, idx(1:step:end),'filled');
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    hca.ZLabel.String = 'v_z (km/s)';
    hs.MarkerFaceAlpha = 0.1;
    hs.MarkerEdgeAlpha = 0.1;
    
    cmap = colormap;
    icolors = round(interp1(1:size(cmap,1),1:size(cmap,1),linspace(1,size(cmap,1),nGroups)));
    group_colors = cmap(icolors,:);
    hca.Title.String = 'k-means';
  end

  if 1 % f(x,y), contour plots of all groups together, gmm
    hca = h(isub); isub = isub + 1;
    holdon = 0;
    for iGroup = 1:numel(pdist_group_gmm)
      %set(hca,'colororder',group_colors(iGroup,:))
      if not(holdon); hold(hca,'on'); holdon = 1; end   
      pd = pdist_group_gmm{iGroup};
      vdf = pd.reduce('2D',[1 0 0],[0 1 0]);
      data = squeeze(vdf.data);
      %data(isnan(data)) = 0;
      nSmooth = 0;
      data = smooth2(data,nSmooth,nSmooth);
      %vdf.data(vdf.data==0) = NaN;
      %vdf.data(vdf.data) = NaN;
      data = log10(data);
  
  
      [hh_,hh] = contour(hca,vdf.depend{1},vdf.depend{1},data',1,'color',group_colors(iGroup,:),'linewidth',linewidth_contour);
      %[hh_,hh] = contour(vdf.depend{1},vdf.depend{1},data',3,'k');
      %hh.FaceColor = 'none';
      %hh.Color = group_colors(iGroup,:);
      %hh.FaceAlpha = 0.2;
      %hh.EdgeColor = group_colors(iGroup,:);

    end
    hold(hca,'off');
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    hold(hca,'on')
    for iGroup = 1:nGroups
      plot(hca,MP_grouped_gmm{iGroup}.sum_vx,MP_grouped_gmm{iGroup}.sum_vy,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
    end
    hold(hca,'off')
  end
  if 1 % f(x,y), contour plots of all groups together, kmeans
    hca = h(isub); isub = isub + 1;
    holdon = 0;
    for iGroup = 1:numel(pdist_group_km)
      %set(hca,'colororder',group_colors(iGroup,:))
      if not(holdon); hold(hca,'on'); holdon = 1; end   
      pd = pdist_group_km{iGroup};
      vdf = pd.reduce('2D',[1 0 0],[0 1 0]);
      data = squeeze(vdf.data);
      %data(isnan(data)) = 0;
      nSmooth = 0;
      data = smooth2(data,nSmooth,nSmooth);
      %vdf.data(vdf.data==0) = NaN;
      %vdf.data(vdf.data) = NaN;
      data = log10(data);
  
  
      [hh_,hh] = contour(hca,vdf.depend{1},vdf.depend{1},data',1,'color',group_colors(iGroup,:),'linewidth',linewidth_contour);
      %[hh_,hh] = contour(vdf.depend{1},vdf.depend{1},data',3,'k');
      %hh.FaceColor = 'none';
      %hh.Color = group_colors(iGroup,:);
      %hh.FaceAlpha = 0.2;
      %hh.EdgeColor = group_colors(iGroup,:);
    end
    hold(hca,'off');
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    hold(hca,'on')    
    for iGroup = 1:nGroups
      plot(hca,MP_grouped_km{iGroup}.sum_vx,MP_grouped_km{iGroup}.sum_vy,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
    end
    hold(hca,'off')
  end

  if 1 % f(x,y), contour plots of all groups together, gmm
    hca = h(isub); isub = isub + 1;
    holdon = 0;
    for iGroup = 1:numel(pdist_group_gmm)
      %set(hca,'colororder',group_colors(iGroup,:))
      if not(holdon); hold(hca,'on'); holdon = 1; end   
      F = gm_out.F{iGroup};
      X = gm_out.X(:,:,1);
      Y = gm_out.Y(:,:,1);
      Fxy = sum(F,3);
      contour(hca,X,Y,Fxy,'color',group_colors(iGroup,:))
    end
    hold(hca,'off');
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
    hold(hca,'on')
    for iGroup = 1:nGroups
      plot(hca,MP_grouped_gmm{iGroup}.sum_vx,MP_grouped_gmm{iGroup}.sum_vy,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
    end
    hold(hca,'off')
  end
  
  
  
  compact_panels(h,0.05)
  c_eval('h(?).XLim = vlim*[-1 1]; h(?).YLim = vlim*[-1 1];',1:numel(h))
  c_eval('h(?).Box = ''on'';',1:numel(h))
  c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))  
  c_eval('h(?).XTick = h(?).YTick;',1:numel(h))
  c_eval('axis(h(?),''square'');',1:numel(h))

  hca = h(isub); isub = isub + 1;
  delete(hca)

  %% Print results
  tint_print = pdist.time.utc('MMHHSS_mmm');
  cn.print(sprintf('gmm_kmeans_%s',tint_print))





end