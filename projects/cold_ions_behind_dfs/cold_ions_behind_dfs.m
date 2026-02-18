
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

db_table_df = db_table_ff(db_table_ff.is_df==1,:);

t0_ff = EpochTT('2017-05-05T19:41:44.790324');
% time: microseconds since 2017-05-05T19:41:44.790324
time_ff_epochtt = irf_time((t0_ff.ttns + db_table_ff.time*1e3),'ttns>EpochTT');



t0_df = EpochTT('2017-05-19T03:06:44.458185978');
% t_df: nanoseconds since  2017-05-19T03:06:44.458185978
time_df_epochtt = irf_time((t0_df.ttns + int64(db_table_df.t_df)),'ttns>EpochTT');


%%
units = irf_units;
% Extract time interval
tint = irf.tint('2017-07-25T22:09:30.00Z/2017-07-25T22:11:00.00Z');
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
c_eval('PD_counts = iPDist?_counts;',ic)
PD_clean = PD_orig.remove_noise(nMean,nThresh,PD_counts);
%PD_diff = PD_orig+-PD_clean;

tsElow = PD_orig.find_noise_energy_limit(5).movmean(30);
emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit

PD_orig_notmasked = PD_orig;
PD_clean_notmasked = PD_clean;
%PD_diff_notmasked = PD_diff;

PD_orig = PD_orig.mask({emask_mat});
PD_clean = PD_clean.mask({emask_mat});
%PD_diff = PD_diff.mask({emask_mat});

PD = PD_orig;
PD = PD_clean;
%% Run through distributions
%time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');

nMovMean = 5; % Number of distributions for the moving average
elim = [200 Inf]; % Energy limit
N_MP = 100000; % Number of macroparticles
nGroupsMax = 2;  % number of classes/groups for kmeans and gmm


% Decide how many output times
dt_all = -5:0.5:5;
%dt_all = 0;
%times = time_xline + dt_all;
dt = 0.5;
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
pdist_all = PD.movmean(nMovMean);
doPrint = 0;

T = dt;
% Run through all times
output = struct();
its = 1:5:pdist_all.length; % 101:102; %100:101%:5:pdist_all.length%1:10:pdist_all.length% 2:(nT-1)%2:(nT-1)
nT = numel(its);
AIC_all = zeros(nT,nGroupsMax);
BIC_all = zeros(nT,nGroupsMax);
    count = 0; % time counter
for it = its
    count = count + 1;
  for nGroups = 1:nGroupsMax  
    %dt = dt_all(it);
    %time = time_xline + dt;  
    %time = times(it);
    %pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).tlim(time + 0.15*0.5*[-1 1]).elim([300 Inf]);
    %pdist = pdist_all.tlim(time + dt*0.5*[-1 1]);
    pdist = pdist_all(it);
    time = pdist.time;
    
    % These should be in dbcs
    %c_eval('density = mean(ne?.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);',ic)
    %c_eval('vel = mean(mvaVi?.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);',ic)
    %c_eval('pressure = mean(mvaPi?.tlim(time + nMovMean* 0.15*0.5*[-1 1]).data,1);',ic)
    c_eval('scpot = mean(scPot?.tlim(time + T* 0.15*0.5*[-1 1]).data,1);',ic)
    scpot = irf.ts_scalar(time,scpot);
   
    % Initial macroparticles, 
    MP = pdist.elim([50 Inf]).macroparticles('ntot',N_MP,'skipzero',1,'scpot',scpot);
    nMP = numel(MP.dv);
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
    
    
    if 0
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
  
    %MP_grouped_vy = macroparticle_moments(MP,idx_vy);
    %pdist_group_vy = partial_pdist(pdist,MP,idx_vy);
  
    %MP_grouped_km = macroparticle_moments(MP,idx_km);
    %pdist_group_km = partial_pdist(pdist,MP,idx_km);
    
  %  MP_grouped_dgb = macroparticle_moments(MP,idx_dgb);
  %  pdist_group_dgb = partial_pdist(pdist,MP,idx_dgb);
  
    % The cuts between kmeans clusters are very (unphyscially) sharp.
    % Refit the separated kmeans clusters to 'recapture' their tails or
    % generally improve the identification. Each kmeans cluster are fitted to
    % a single normal distribution.
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
    
    
    gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],nGroups,'Start',S,'SharedCovariance',false);
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
  
    AIC_all(count,nGroups) = gm_out.gm.AIC;
    BIC_all(count,nGroups) = gm_out.gm.BIC;
  
    %MP_grouped_gmm = macroparticle_moments(MP,gm_out.idx);
    %pdist_group_gmm = partial_pdist(pdist,MP,gm_out.idx);
    
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
  
    % Get the tempreature of each component 
    vt2 = gm_out.gm.Sigma*1e6; % (km/s)^2 -> (m/s)^2
    cov_T_mat = units.mp*vt2/2/units.eV;
    clear T_tmp
    for iGroup = 1:nGroups
      T_tmp(iGroup) = trace(cov_T_mat(:,:,iGroup))/3;
    end
  
    % Get the density of each component
    n_scale = sum(MP.df.*MP.dv); % cc
  
    % Get the velocity of each component
    clear vx_tmp  vy_tmp vz_tmp
    for iGroup = 1:nGroups
      vx_tmp(iGroup) = gm_out.gm.mu(iGroup,1);
      vy_tmp(iGroup) = gm_out.gm.mu(iGroup,2);
      vz_tmp(iGroup) = gm_out.gm.mu(iGroup,3);
    end
  
  
  
    % Gather output
  
    % Posterior probabilite histogram/pdf
    clear N
    edges_post_prob = 0:0.02:1;
    for iGroup = 1:nGroups    
      N_tmp = histcounts(P(:,iGroup),edges_post_prob,'Normalization','probability');
      output(count,nGroups).PosteriorPDF{iGroup} = N_tmp;
    end
  
    output(count,nGroups).time = time;
    output(count,nGroups).n = gm_out.gm.ComponentProportion*n_scale;
    output(count,nGroups).vx = vx_tmp;
    output(count,nGroups).vy = vy_tmp;
    output(count,nGroups).vz = vz_tmp;
    output(count,nGroups).T = T_tmp;
    output(count,nGroups).nlogT = gm_out.gm.NegativeLogLikelihood;
    output(count,nGroups).AIC = gm_out.gm.AIC;
    output(count,nGroups).BIC = gm_out.gm.BIC; 
  
  
    if 0 % Plot
      %%
    vdf_12 = pdist.reduce('2D',[1 0 0],[0 1 0]);
    vdf_1 = pdist.reduce('1D',[1 0 0]);
  
    h = setup_subplots(3,2);
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
    %data = log10(sum(gm_out.F{2},3));
    data = sum(gm_out.F{1},3);
    for iGroup = 2:nGroups
      data = data + sum(gm_out.F{iGroup},3);
    end
    pcolor(hca,gm_out.vx,gm_out.vy,log10(data)')
    shading(hca,'flat')
    
    hca = h(isub); isub = isub + 1;
    clear data
    for iGroup = 1:nGroups
      data(iGroup,:) = sum(gm_out.F{iGroup},[2 3]);
    end
    % Scale to density: sum(MP.df.*MP.dv)*1e-6
    dv = (gm_out.vx(2)-gm_out.vx(1))*1e3;
    plot(hca,vdf_1.depend{1},vdf_1.data,gm_out.vx,sum(data,1)*n_scale*1e-6*dv*dv*1e3)
    shading(hca,'flat')
  
    hca = h(isub); isub = isub + 1;
    histogram(hca,P(:,1),0:0.02:1,'Normalization','probability')
    legs_hist = {'G1'};
    hold(hca,'on')
    for iGroup = 2:nGroups-0
      histogram(hca,P(:,iGroup),0:0.02:1,'Normalization','probability')
      legs_hist{iGroup} = sprintf('G%g',iGroup);
    end
    hold(hca,'off')
    hca.YLim = [0 0.05];
    hca.XLabel.String = 'Posterior probability';
    hca.YLabel.String = 'Probability';
    legend(hca,legs_hist,'location','north')
    
  
  
    c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on''; h(?).Layer = ''top'';',1:numel(h))
    linkprop(h(1:4),{'XLim','YLim'})
    %compact_panels(h,0.02,0.1)
    pause(0.2)
    end  
  end
  %output(:,1).AIC
end

%% Plot results as timeseries
times = cat(1,[output(:,1).time]);
%for it = 1:nT
%  for iGroup = 1:nGroupsMax
%    data_T(it,iGroup) = output(it,iGroup).T;
%  end
%end
iGroup = 2;
tsT = irf.ts_scalar(times,cat(1,output(:,iGroup).T));
tsN = irf.ts_scalar(times,cat(1,output(:,iGroup).n));
tsVx = irf.ts_scalar(times,cat(1,output(:,iGroup).vx));
tsVy = irf.ts_scalar(times,cat(1,output(:,iGroup).vy));
tsVz = irf.ts_scalar(times,cat(1,output(:,iGroup).vz));
tsNlogT = irf.ts_scalar(times,cat(1,output(:,iGroup).nlogT));
tsAIC = irf.ts_scalar(times,cat(1,output(:,iGroup).AIC));
tsBIC = irf.ts_scalar(times,cat(1,output(:,iGroup).BIC));

if 0
%clear NP
%for iGroup = 1:nGroups
%  for it = 1:numel(output)
%    NP(it,:,iGroup) = output(it,iGroup).PosteriorPDF{iGroup};
%  end
%end

clear specrec specrec_tmp
for iGroup = 1:nGroups  
  d_edges = diff(edges_post_prob);
  specrec_tmp.f = edges_post_prob(2:end)-d_edges*0.5;
  specrec_tmp.p = NP(:,:,iGroup);
  specrec_tmp.t = times.epochUnix;
  specrec_tmp.f_label = {'Posterior','probability'};
  specrec_tmp.p_label = {'Probability'};
  specrec{iGroup} = specrec_tmp;
end
end

h = irf_plot(9);


if 1
  hca = irf_panel('AIC');
  irf_plot(hca,tsAIC,'*')
  hca.YLabel.String = {'AIC'};
  
  hca = irf_panel('BIC');
  irf_plot(hca,tsBIC,'*')
  hca.YLabel.String = {'BIC'};
  
  hca = irf_panel('nlogT');
  irf_plot(hca,tsNlogT,'*')
  hca.YLabel.String = {'Negative','Log Likelihood'};  
end

if 0 % posterior probability
  for iGroup = 1:nGroups  
    hca = irf_panel(sprintf('posterior probability G%g',iGroup));
    %hca.YLabel.String = 'n (cc)';
    irf_spectrogram(hca,specrec{iGroup},'log')    
  end
end

hca = irf_panel('n');
irf_plot(hca,tsN,'*')
hca.YLabel.String = 'n (cc)';

hca = irf_panel('vx');
irf_plot(hca,tsVx,'*')
hca.YLabel.String = 'v_x (km/s)';

hca = irf_panel('vx');
irf_plot(hca,tsVx,'*')
hca.YLabel.String = 'v_x (km/s)';

hca = irf_panel('vy');
irf_plot(hca,tsVy,'*')
hca.YLabel.String = 'v_y (km/s)';

hca = irf_panel('vz');
irf_plot(hca,tsVz,'*')
hca.YLabel.String = 'v_z (km/s)';

hca = irf_panel('T');
irf_plot(hca,tsT,'*')
hca.YLabel.String = 'T (eV)';

hca = irf_panel('ion deflux omni');
%hca.YLabel.String = 'n (cc)';
irf_spectrogram(hca,PD.deflux.omni.specrec,'log')
hca.YScale = "log";

irf_plot_axis_align
irf_zoom(h,'x',[PD.time.start PD.time.stop])
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))


%% Compare with thu enu



hca = subplot(4,1,1);
nLogT_all=[cat(1,output(:,1).nlogT),cat(1,output(:,2).nlogT),cat(1,output(:,3).nlogT),cat(1,output(:,4).nlogT)];
plot(hca,nLogT_all);

hca = subplot(4,1,2);
plot(hca,AIC_all/nMP);

hca = subplot(4,1,3);
plot(hca,BIC_all/nMP);

hca = subplot(4,1,4);
plot(hca,BIC_all./repmat(BIC_all(:,1),[1 nGroupsMax]));