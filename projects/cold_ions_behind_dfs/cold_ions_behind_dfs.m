localuser = 'cecilianorgren';
file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_b_gsm_2017-2021.nc'];

file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.nc'];
%data = load(file);

%ncdisp(file)

info = ncinfo(file);
vars = {info.Variables.Name};
nvars = numel(vars);

clear db
for ivar = 1:nvars
  db.(vars{ivar}) = ncread(file,vars{ivar});
end

db_table_ff = struct2table(db);

db_table_df = db_table(db_table_ff.is_df==1,:);

t0_ff = EpochTT('2017-05-05T19:41:44.790324');
% time: microseconds since 2017-05-05T19:41:44.790324
time_ff_epochtt = irf_time((t0_ff.ttns + db_table_ff.time*1e3),'ttns>EpochTT');



t0_df = EpochTT('2017-05-19T03:06:44.458185978');
% t_df: nanoseconds since  2017-05-19T03:06:44.458185978
time_df_epochtt = irf_time((t0_df.ttns + int64(db_table_df.t_df)),'ttns>EpochTT');


%%

% Extract time interval
tint = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');

% Load data
ic = 1;
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)


%% Clean PDist
nMean = [5,3,3,3]; nThresh = 3;

c_eval('PD_orig = iPDist?;',ic)
%c_eval('PD_counts = iPDist?_counts;',ic)
%PD_clean = PD_orig.remove_noise(nMean,nThresh,PD_counts);
%PD_diff = PD_orig+-PD_clean;

tsElow = PD_orig.find_noise_energy_limit(5).movmean(30);
emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit

%PD_orig_notmasked = PD_orig;
%PD_clean_notmasked = PD_clean;
%PD_diff_notmasked = PD_diff;

PD_orig = PD_orig.mask({emask_mat});
%PD_clean = PD_clean.mask({emask_mat});
%PD_diff = PD_diff.mask({emask_mat});

PD = PD_orig;
%% Run through distributions
%time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');

nMovMean = 7; % Number of distributions for the moving average
elim = [200 Inf]; % Energy limit
N_MP = 30000; % Number of macroparticles
nGroups = 2;  % number of classes/groups for kmeans and gmm


% Decide how many output times
dt_all = -5:0.5:5;
%dt_all = 0;
%times = time_xline + dt_all;
dt = 1;
times = tint.start:dt:tint.stop;
nT = times.length;

% Initialize arrays for output
for iGroup = 1:nGroups
  out.vy.n{iGroup} = zeros(nT,1);   % density gradient based scan
  out.vy.v{iGroup} = zeros(nT,3);
  out.vy.p{iGroup} = zeros(nT,3,3);
  out.vy.n_MP{iGroup} = zeros(nT,1);   % density gradien based scan, moments from macroparticles
  out.vy.v_MP{iGroup} = zeros(nT,3);
  out.vy.p_MP{iGroup} = zeros(nT,3,3);
  
  out.dgb.n{iGroup} = zeros(nT,1);   % density gradient based scan
  out.dgb.v{iGroup} = zeros(nT,3);
  out.dgb.p{iGroup} = zeros(nT,3,3);
  out.dgb.n_MP{iGroup} = zeros(nT,1);   % density gradien based scan, moments from macroparticles
  out.dgb.v_MP{iGroup} = zeros(nT,3);
  out.dgb.p_MP{iGroup} = zeros(nT,3,3);

  out.km.n{iGroup} = zeros(nT,1);   % k-means
  out.km.v{iGroup} = zeros(nT,3);
  out.km.p{iGroup} = zeros(nT,3,3);
  out.km.n_MP{iGroup} = zeros(nT,1);   % k-means, moments from macroparticles
  out.km.v_MP{iGroup} = zeros(nT,3);
  out.km.p_MP{iGroup} = zeros(nT,3,3);
  
  out.gm.pdf = cell(nT,1);
  out.gm.model = cell(nT,1);
  out.gm.n{iGroup} = zeros(nT,1);   % gaussian mixture model
  out.gm.v{iGroup} = zeros(nT,3);
  out.gm.p{iGroup} = zeros(nT,3,3);
  out.gm.n_MP{iGroup} = zeros(nT,1);   % gaussian mixture model, moments from macroparticles
  out.gm.v_MP{iGroup} = zeros(nT,3);
  out.gm.p_MP{iGroup} = zeros(nT,3,3);

  out.pr.n{iGroup} = zeros(nT,1);   % pressure rotation
  out.pr.v{iGroup} = zeros(nT,3);
  out.pr.p{iGroup} = zeros(nT,3,3);
end

% Prepare the distribution
pdist_all = PD;
doPrint = 0;

T = dt;
% Run through all times
for it = 2%nT-1%2:(nT-1)
  %dt = dt_all(it);
  %time = time_xline + dt;  
  time = times(it);
  %pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).tlim(time + 0.15*0.5*[-1 1]).elim([300 Inf]);
  pdist = pdist_all.tlim(time + dt*0.15*0.5*[-1 1]);
  vdf_12 = pdist.reduce('2D',[1 0 0],[0 1 0]);
  
  % These should be in dbcs
  %c_eval('density = mean(ne?.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);',ic)
  %c_eval('vel = mean(mvaVi?.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);',ic)
  %c_eval('pressure = mean(mvaPi?.tlim(time + nMovMean* 0.15*0.5*[-1 1]).data,1);',ic)
  c_eval('scpot = mean(scPot?.tlim(time + T* 0.15*0.5*[-1 1]).data,1);',ic)
  scpot = irf.ts_scalar(time,scpot);
 
  % Initial macroparticles, 
  MP = pdist.elim([50 Inf]).macroparticles('ntot',N_MP,'skipzero',1,'scpot',scpot);
  %nfrac = sum(MP.dv.*MP.df)*1e-6/mean(ne3.tlim(time + 0.15*0.5*[-1 1]).data,1);
  MP.dn = MP.df.*MP.dv;%/nfrac;
  V_dbcs = [MP.vx, MP.vy, MP.vz]; 
  
  % Need to rotate these into the specified coordinate system  
  % The rotation is not done exactly right, due to not taking into account
  % differences between dbcs and gse. To do in the future
  %V = V_dbcs*lmn';
  V = V_dbcs;
  MP.vx = V(:,1);
  MP.vy = V(:,2);
  MP.vz = V(:,3);

  % Initial guess for 4 populations
  initial_guess = [0 -1000 -1000;...
                   0 -1000 +1000;...
                   0 +1000 -1000;...
                   0 +1000 +1000;...
                   0  0000  0000;...
                   0  0000  0000];
  initial_guess = initial_guess(1:nGroups,:);
  
  % kmeans  
  [idx_km, c, sumd, d] = kmeans(V, nGroups,'Start',initial_guess);
  %s = silhouette(V, idx, 'sqeuclid');

  % dgb
  %cluster1 = clusterDBSCAN('MinNumPoints',N_MP/10,'Epsilon',2, ...
  %  'EnableDisambiguation',false);
  %[idx_dgb,clusterids] = cluster1(X)
  
  
  if 1
  [ccount eedg mmid lloc] = histcn([MP.vy,MP.vz],[-Inf 0 Inf],[-Inf 0 Inf]);
  
  idx_vy_ = zeros(numel(MP.dn),4);
  idx_vy_(:,1) = (lloc(:,1)==1).*(lloc(:,2)==1);
  idx_vy_(:,2) = (lloc(:,1)==1).*(lloc(:,2)==2);
  idx_vy_(:,3) = (lloc(:,1)==2).*(lloc(:,2)==1);
  idx_vy_(:,4) = (lloc(:,1)==2).*(lloc(:,2)==2);
  
  idx_vy = zeros(numel(MP.dn),1);
  idx_vy(find(idx_vy_(:,1)==1)) = 1;
  idx_vy(find(idx_vy_(:,2)==1)) = 2;
  idx_vy(find(idx_vy_(:,3)==1)) = 3;
  idx_vy(find(idx_vy_(:,4)==1)) = 4;

  else 
    idx_negvy = find(MP.vy<0);
    idx_posvy = find(MP.vy>0);
    idx_vy = zeros(numel(MP.dn),1);
    idx_vy(idx_negvy) = 1;
    idx_vy(idx_posvy) = 2;
  end

  MP_grouped_vy = macroparticle_moments(MP,idx_vy);
  pdist_group_vy = partial_pdist(pdist,MP,idx_vy);

  MP_grouped_km = macroparticle_moments(MP,idx_km);
  pdist_group_km = partial_pdist(pdist,MP,idx_km);
  
%  MP_grouped_dgb = macroparticle_moments(MP,idx_dgb);
%  pdist_group_dgb = partial_pdist(pdist,MP,idx_dgb);

  % The cuts between kmeans clusters are very (unphyscially) sharp.
  % Refit the separated kmeans clusters to 'recapture' their tails or
  % generally improve the identification. Each kmeans cluster are fitted to
  % a single normal distribution.
  %
  if 0
  for iGroup = 1:nGroups
    S.mu = [MP_grouped_km{iGroup}.sum_vx, MP_grouped_km{iGroup}.sum_vy, MP_grouped_km{iGroup}.sum_vz];
    S.Sigma = [200 10 10; 10 200 10; 10 10 200]; % these numbers are a bit arbitrary
    S.ComponentProportion = 1;      
        
    Vtmp = V(idx_km==iGroup,:);
    gm_out_kmeans{iGroup} = gaussian_mixture_model(Vtmp,1,'Start',S);    
    
    
    if 0 % plot results   
      %%
      % Setup grid
      vmax = 2400;
      nv = 80;
      vvec = linspace(-vmax,vmax,nv);      
      dv = (vvec(2)-vvec(1))*1e3;
      [X,Y,Z] = ndgrid(vvec,vvec,vvec);
      XYZ = [X(:) Y(:) Z(:)];
      
      % Extract parameters
      mu = gm_out_kmeans{iGroup}.gm.mu;
      Sigma = gm_out_kmeans{iGroup}.gm.Sigma;
      Ftmp = gm_out_kmeans{iGroup}.gm.ComponentProportion*mvnpdf(XYZ, mu, Sigma); 
      F = reshape(Ftmp,size(X));
      vdf_model_x = squeeze(sum(sum(F,3),2)); % 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,3),2));
      vdf_model_y = squeeze(sum(sum(F,3),1)); % 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,3),1));
      vdf_model_z = squeeze(sum(sum(F,1),2)); % 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,1),2));
    
      histvec = -vmax:100:vmax;
      h = setup_subplots(3,1); % nDim x 1
      isub = 1;
      
      hca = h(isub); isub = isub + 1;
      histogram(hca,Vtmp(:,1),histvec)
      hold(hca,'on')
      yyaxis(hca,'right')
      plot(hca,vvec,vdf_model_x,'k')
      hold(hca,'off')
      hca.XLabel.String = 'v_L (km/s)';
      
      
      hca = h(isub); isub = isub + 1;
      histogram(hca,Vtmp(:,2),histvec)
      hold(hca,'on')
      yyaxis(hca,'right')
      plot(hca,vvec,vdf_model_y,'k')
      hold(hca,'off')
      hca.XLabel.String = 'v_M (km/s)';
      
      
      hca = h(isub); isub = isub + 1;
      histogram(hca,Vtmp(:,3),histvec)
      hold(hca,'on')
      yyaxis(hca,'right')
      plot(hca,vvec,vdf_model_z,'k')
      hold(hca,'off')
      hca.XLabel.String = 'v_N (km/s)';
      
    end
        
  end
  end
  if 0 % plot results comparing the four clusters
     %%
      % Setup grid
      vmax = 2400;
      nv = 100;
      vvec = linspace(-vmax,vmax,nv);      
      dv = (vvec(2)-vvec(1))*1e3;
      [X,Y,Z] = ndgrid(vvec,vvec,vvec);
      XYZ = [X(:) Y(:) Z(:)];
      
      % Extract parameters
      mu = gm_out_kmeans{iGroup}.gm.mu;
      Sigma = gm_out_kmeans{iGroup}.gm.Sigma;
      Ftmp = gm_out_kmeans{iGroup}.gm.ComponentProportion*mvnpdf(XYZ, mu, Sigma); 
      F = reshape(Ftmp,size(X));
      vdf_model_x = squeeze(sum(sum(F,3),2)); % 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,3),2));
      vdf_model_y = squeeze(sum(sum(F,3),1)); % 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,3),1));
      vdf_model_z = squeeze(sum(sum(F,1),2)); % 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,1),2));
    
      histvec = -vmax:50:vmax;
      % Bin particles
      [count, edges, mid, loc] = histcn([V idx_km],histvec,histvec,histvec,0.5:1:4.5);
      
      h = setup_subplots(3,1); % nDim x 1
      isub = 1;
      
      if 1 % vL
        hca = h(isub); isub = isub + 1;            
        count_v = squeeze(sum(sum(count(:,:,:,:),3),2));
        plot(hca,mid{1},sum(count_v,2),'k');
        hold(hca,'on')
        plot(hca,mid{1},count_v);
        %for iGroup = 1:4
          %hbar = bar(hca,mid{1},count_vx(:,iGroup),1,'facealpha',0.5);
        %end
        hold(hca,'off')      
        hca.XLabel.String = 'v_L (km/s)';
      end
      if 1 % VM
        hca = h(isub); isub = isub + 1;            
        count_v = squeeze(sum(sum(count(:,:,:,:),3),1));
        plot(hca,mid{1},sum(count_v,2),'k');
        hold(hca,'on')
        plot(hca,mid{1},count_v);
        %for iGroup = 1:4
          %hbar = bar(hca,mid{1},count_vx(:,iGroup),1,'facealpha',0.5);
        %end
        hold(hca,'off')      
        hca.XLabel.String = 'v_M (km/s)';
      end
      if 1 % VN
        hca = h(isub); isub = isub + 1;            
        count_v = squeeze(sum(sum(count(:,:,:,:),2),1));
        plot(hca,mid{1},sum(count_v,2),'k');
        hold(hca,'on')
        plot(hca,mid{1},count_v);
        %for iGroup = 1:4
          %hbar = bar(hca,mid{1},count_vx(:,iGroup),1,'facealpha',0.5);
        %end
        hold(hca,'off')      
        hca.XLabel.String = 'v_N (km/s)';
      end
      %%
      hca = h(isub); isub = isub + 1;
      histogram(hca,Vtmp(:,2),histvec)
      hold(hca,'on')
      yyaxis(hca,'right')
      plot(hca,vvec,vdf_model_y,'k')
      hold(hca,'off')
      hca.XLabel.String = 'v_M (km/s)';
      
      
      hca = h(isub); isub = isub + 1;
      histogram(hca,Vtmp(:,3),histvec)
      hold(hca,'on')
      yyaxis(hca,'right')
      plot(hca,vvec,vdf_model_z,'k')
      hold(hca,'off')
      hca.XLabel.String = 'v_N (km/s)';
  
  end
  %
  % Gaussian mixture model
  % Use results from kmeans to give initial guess of gmm
  
  initial_guess_method = 'guess1';
  switch initial_guess_method
    case 'kmeans'
      for iGroup = 1:nGroups
        S.mu(iGroup,:) = [MP_grouped_km{iGroup}.sum_vx, MP_grouped_km{iGroup}.sum_vy, MP_grouped_km{iGroup}.sum_vz];
        S.Sigma(1:3,1:3,iGroup) = [200 10 10; 10 200 10; 10 10 200];
        S.ComponentProportion(iGroup) = 1;                 
      end
    case 'guess1'  
      S.mu = initial_guess; % same as for kmeans         
      S.Sigma = repmat([500 10 10; 10 500 10; 10 10 500],[1 1 nGroups]);  
      S.ComponentProportion = repmat(1,[1,nGroups]);
    case 'guess2' % a bit colder distributions 
      S.mu = [0 -1000 -1000;...
              0 -1000 +1000;...
              0 +100 -1000;...
              0 +100 +1000;...
              0 0 0;...
              0 0 0];
      S.mu = S.mu(1:nGroups,:);
      S.Sigma = repmat([100 10 10; 10 100 10; 10 10 100],[1 1 nGroups]);  
      S.ComponentProportion = repmat(1,[1,nGroups]);
  end
  
  
  gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],nGroups,'Start',S);
  %gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],nGroups,'Replicates',5);
  %gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],nGroups,'Start',S,'CovType','diagonal');
  %gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],nGroups);
  % Estimate goodness of fit
  goodness.mahal = mahal(gm_out.gm, gm_out.R);
  goodness.posterior = posterior(gm_out.gm, gm_out.R);
  [P,nlogL] = posterior(gm_out.gm, gm_out.R);
  goodness.NegativeLogLikelihood = gm_out.gm.NegativeLogLikelihood;
  goodness.AIC = gm_out.gm.AIC;
  goodness.BIC = gm_out.gm.BIC;

  MP_grouped_gmm = macroparticle_moments(MP,gm_out.idx);
  pdist_group_gmm = partial_pdist(pdist,MP,gm_out.idx);
  
  % Pressure-based rotation
  % Rotate distribution in the xy-plane so that the pressure is maximized
  % along one dimension. Then 
%  [T,vdf] = rotate_vdf_based_on_pressure(pdist,pressure,vel);
  % Fit 1D distributions
  % ...
  
  % Collect output
  if 0
  for iGroup = 1:nGroups
    if 1%iGroup<3
      out.vy.n_MP{iGroup}(it,1) = MP_grouped_vy{iGroup}.sum_n;
      out.vy.v_MP{iGroup}(it,1:3) = [MP_grouped_vy{iGroup}.sum_vx, MP_grouped_vy{iGroup}.sum_vy, MP_grouped_vy{iGroup}.sum_vz];
    end
    
    % k-means clustering
    out.km.n_MP{iGroup}(it,1) = MP_grouped_km{iGroup}.sum_n;
    out.km.v_MP{iGroup}(it,1:3) = [MP_grouped_km{iGroup}.sum_vx, MP_grouped_km{iGroup}.sum_vy, MP_grouped_km{iGroup}.sum_vz];
    
    
    % Gaussian mixture model
    out.gm.pdf{it} = gm_out.gmPDF;
    out.gm.model{it} = gm_out.gm;
    out.gm.n{iGroup}(it,1) = gm_out.gm.ComponentProportion(iGroup)*density; % needs to be multiplied with total density
    out.gm.v{iGroup}(it,1:3) = gm_out.gm.mu(iGroup,:);
    out.gm.p{iGroup}(it,:,:) = gm_out.gm.Sigma(:,:,iGroup); % for now, only returns diagonal Sigma (pressure), off-diag implemented
    %out.gm.p{iGroup}(it,2) = gm_out.gm.Sigma(1,2,iGroup);
    %out.gm.p{iGroup}(it,3) = gm_out.gm.Sigma(1,3,iGroup);

    out.gm.n_MP{iGroup}(it,1) = MP_grouped_gmm{iGroup}.sum_n;
    out.gm.v_MP{iGroup}(it,1:3) = [MP_grouped_gmm{iGroup}.sum_vx, MP_grouped_gmm{iGroup}.sum_vy, MP_grouped_gmm{iGroup}.sum_vz];
    %out.gm.p_MP{iGroup}(it,1) = MP_grouped_gmm{iGroup}.sum_n;
  end
  end

  % Plot comparison
  %fig = figure(71);
  %h = setup_subplots(1,1);
  
  if 0 % plot f(vy) for different pops
    [counts edges mid loc] = histcn([MP.vy,idx_km],[-2000:50:2000],[0.5:1:4.5],'accumdata',MP.dn);
    plot(mid{1},counts,mid{1},sum(counts(:,1:2),2),'--',mid{1},sum(counts(:,3:4),2),'--',mid{1},sum(counts(:,1:4),2),'k:','linewidth',2)
    xlabel('v_y')
    legend({'1','2','3','4','1-2','3-4','1-4'},'location','northeast')
    title(time.utc('yyyy-dd-mm HH:MM:SS.mmm'))
    cn.print(['ion_dyn_kmeans_pops_1234_' time.utc('yyyy-dd-mm_HHMMSS_mmm')])
    pause(0.01)
  end


  h = setup_subplots(2,2);
  isub = 1;
  
  hca = h(isub); isub = isub + 1;
  vdf_12.plot_plane(hca);

  hca = h(isub); isub = isub + 1;
  scatter(hca,MP.vx,MP.vy,5,gm_out.idx)

  hca = h(isub); isub = isub + 1;
  data = log10(sum(gm_out.F{1},3));
  pcolor(hca,gm_out.vx,gm_out.vy,data')
  shading(hca,'flat')

  hca = h(isub); isub = isub + 1;
  data = log10(sum(gm_out.F{2},3));
  data = log10(sum(gm_out.F{1},3)+sum(gm_out.F{2},3));
  pcolor(hca,gm_out.vx,gm_out.vy,data')
  shading(hca,'flat')
  
  c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Layer = ''top'';',1:numel(h))
  linkprop(h,{'XLim','YLim'})
  compact_panels(h,0.01,0.01)
  pause(0.2)
end

% Make timeseries of results
c_eval('kmN?_MP = irf.ts_scalar(times,out.km.n_MP{?});',1:nGroups)
kmN1234_MP_ = kmN1_MP + kmN2_MP + kmN3_MP + kmN4_MP;
c_eval('kmV?_MP = irf.ts_vec_xyz(times,out.km.v_MP{?});',1:nGroups)



c_eval('vyN?_MP = irf.ts_scalar(times,out.vy.n_MP{?});',1:nGroups)
%kmN1234_MP_ = kmN1_MP + kmN2_MP + kmN3_MP + kmN4_MP;
c_eval('vyV?_MP = irf.ts_vec_xyz(times,out.vy.v_MP{?});',1:nGroups)
vyN1234_MP_ = vyN1_MP + vyN2_MP + vyN3_MP + vyN4_MP;
vyV_MP = (vyN1_MP*vyV1_MP + vyN2_MP*vyV2_MP + vyN3_MP*vyV3_MP + vyN4_MP*vyV4_MP)/vyN1234_MP_;

c_eval('gmN? = irf.ts_scalar(times,out.gm.n{?});',1:nGroups)
c_eval('gmV? = irf.ts_vec_xyz(times,out.gm.v{?});',1:nGroups)
c_eval('gmP? = irf.ts_tensor_xyz(times,out.gm.p{?});',1:nGroups)

c_eval('gmN?_MP = irf.ts_scalar(times,out.gm.n_MP{?});',1:nGroups)
c_eval('gmV?_MP = irf.ts_vec_xyz(times,out.gm.v_MP{?});',1:nGroups)


%gmN = irf.ts_scalar(times,out.gm.n);
%gmV = irf.ts_vec_xyz(times,out.gm.v);



