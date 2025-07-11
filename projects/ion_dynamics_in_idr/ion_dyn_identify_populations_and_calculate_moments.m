% Identify different ion populations and their moments for the Torbert
% event.

%% Load data
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
mms.db_init('local_file_db','/Volumes/mms');
%db_info = datastore('mms_db');

units = irf_units;

ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); % Torbert event 

c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
%c_eval('gseB? = mms.get_data(''B_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

%% Calculate moments from reduced vdf
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');
tint_red = time_xline + [-20 20];
%c_eval('if?y = iPDist?.tlim(tint_red).elim([200 Inf]).reduce(''1D'',[0 1 0]);',ic)

vv = if3y.depend{1}(1,:);
dv = vv(2) - vv(1);
nv = numel(vv);
vv_edge = [vv(1:end)-0.5*dv vv(end) + 0.5*dv];
n1 = sum(if3y.data(:,1:nv/2)*dv*1e3,2)*1e-6;
n2 = sum(if3y.data(:,(nv/2+1):nv)*dv*1e3,2)*1e-6;

%VV = repmat(vv,if3y.length,1);
%JJ = if3y.data.*VV*1e3;
%j1 = sum(JJ(:,1:nv/2),2)*dv*1e3*1e-6;
%j2 = sum(JJ(:,(nv/2+1):nv),2)*dv*1e3*1e-6;



tsN1 = irf.ts_scalar(if3y.time,n1);
tsN2 = irf.ts_scalar(if3y.time,n2);
tsN = ne3.tlim(tint_red).resample(if3y);


%for it = 1%:if3y.length
  c_eval('pdist = iPDist?.tlim(tint_red).elim([200 Inf]);',ic)
  %vdf = pdist.tlim(tint_red).elim([200 Inf]).reduce('2D',[1 0 0],[0 1 0]);
  vv = vdf.depend{1}(1,:);
  dv = vv(2) - vv(1);
  nv = numel(vv);
  vv_edge = [vv(1:end)-0.5*dv vv(end) + 0.5*dv];  
  [VVx,VVy] = ndgrid(1:vdf.length,vv,vv);  
  [VVx_,VVy] = ndgrid(vv,vv);  
  VVx_ = repmat(VVx_,1,1,vdf.length);
  VVx = permute(VVx_,[3 1 2]);

  NNx = vdf.data.*(dv*1e3)^2;  
  NNx1 = sum(NNx(:,:,1:nv/2),3);
  NNx2 = sum(NNx(:,:,(nv/2+1):nv),3);
  n1 = sum(NNx1,2);
  n2 = sum(NNx2,2);

  JJx = vdf.data.*VVx*1e3*(dv*1e3)^2; 
  JJx1 = sum(JJx(:,:,1:nv/2),3);
  JJx2 = sum(JJx(:,:,(nv/2+1):nv),3);
  j1 = sum(JJx1,2);
  j2 = sum(JJx2,2);
  
  tsN1 = irf.ts_scalar(vdf.time,n1*1e-6);
  tsN2 = irf.ts_scalar(vdf.time,n2*1e-6);
  tsN = ne3.tlim(tint_red).resample(if3y);

  
  tsJ1 = irf.ts_scalar(vdf.time,j1);
  tsJ2 = irf.ts_scalar(vdf.time,j2);
  
  tsV1 = irf.ts_scalar(vdf.time,j1./n1*1e-3);
  tsV2 = irf.ts_scalar(vdf.time,j2./n2*1e-3);

  tsV1frac = irf.ts_scalar(vdf.time,n1.*j1./(n1+n2));
  tsV2frac = irf.ts_scalar(vdf.time,n2.*j2./(n1+n2));
  
  if 1
    h = irf_plot(4);

    hca = irf_panel('n');
    hca.ColorOrder = pic_colors('matlab');
    irf_plot(hca,{tsN1,tsN2,tsN1+tsN2,ni3.tlim([tsN1.time.start tsN1.time.stop]),ne3.tlim([tsN1.time.start tsN1.time.stop])},'comp')
    irf_legend(hca,{'vy>0','vy<0','all vy','ni fpi','ne fpi'},[0.02 0.98])
    hca.YLabel.String = 'n (cc)';

    hca = irf_panel('j');
    hca.ColorOrder = pic_colors('matlab');
    irf_plot(hca,{tsJ1,tsJ2},'comp')
    irf_legend(hca,{'vy>0','vy<0'},[0.02 0.98])
    hca.YLabel.String = 'j';

    hca = irf_panel('vfrac');
    hca.ColorOrder = pic_colors('matlab');
    irf_plot(hca,{tsV1frac,tsV2frac},'comp')
    irf_legend(hca,{'vy>0','vy<0'},[0.02 0.98])
    hca.YLabel.String = 'n_sv_s/n_{tot}';

    hca = irf_panel('v');
    hca.ColorOrder = pic_colors('matlab');
    irf_plot(hca,{tsV1,tsV2,gseVi3.x},'comp')
    irf_legend(hca,{'vy>0','vy<0','v_i fpi'},[0.02 0.98])
    hca.YLabel.String = 'v = j_s/n_s';

    irf_zoom(h,'x',tsN1.time)
  end

%end

%% Local coordinate system
% From Torbert
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn_edr = [L;M;N];

% Don't remember, vi MVA?
L_vi = -[-0.8906    0.4548    0.0045];
M_vi = [ 0.4539    0.8893   -0.0559];
N_vi = -[-0.0294   -0.0477   -0.9984];
lmn_vi = [L_vi; M_vi; N_vi];

% GSE
L_gse = [1 0 0];
M_gse = [0 1 0];
N_gse = [0 0 1];
lmn_gse = [L_gse; M_gse; N_gse];

% Ion VDF by eye
L_gse = [1 0 -.2]; L_gse = L_gse/norm(L_gse);
M_gse = [0 1 -0.2]; M_gse = cross(L_gse,cross(M_gse,L));
N_gse = cross(L_gse,M_gse);
lmn_fi = [L_gse; M_gse; N_gse];

lmn = lmn_edr;
lmn = lmn_fi;
lmn = lmn_gse;
lmn = lmn_fi;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);

c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
%c_eval('mvaPi? = lmn*gseVi?*lmn''; mvaVi?.name = ''Pi LMN'';',ic) % check rotation
c_eval('mvaPi? = gsePi?''; mvaPi?.name = ''Pi LMN'';',ic) % check rotation


%% Run through distributions
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');

nMovMean = 7; % Number of distributions for the moving average
elim = [200 Inf]; % Energy limit
N_MP = 30000; % Number of macroparticles
nGroups = 4;  % number of classes/groups for kmeans and gmm


% Decide how many output times
dt_all = -5:0.5:5;
%dt_all = 0;
times = time_xline + dt_all;
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
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)
doPrint = 0;

% Run through all times
for it = 1:nT
  %dt = dt_all(it);
  %time = time_xline + dt;  
  time = times(it);
  %pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).tlim(time + 0.15*0.5*[-1 1]).elim([300 Inf]);
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  
  % These should be in dbcs
  c_eval('density = mean(ne?.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);',ic)
  c_eval('vel = mean(mvaVi?.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);',ic)
  c_eval('pressure = mean(mvaPi?.tlim(time + nMovMean* 0.15*0.5*[-1 1]).data,1);',ic)
  c_eval('scpot = mean(scPot?.tlim(time + nMovMean* 0.15*0.5*[-1 1]).data,1);',ic)
  scpot = irf.ts_scalar(time,scpot);
 
  % Initial macroparticles, 
  MP = pdist.elim([50 Inf]).macroparticles('ntot',N_MP,'skipzero',1,'scpot',scpot);
  %nfrac = sum(MP.dv.*MP.df)*1e-6/mean(ne3.tlim(time + 0.15*0.5*[-1 1]).data,1);
  MP.dn = MP.df.*MP.dv;%/nfrac;
  V_dbcs = [MP.vx, MP.vy, MP.vz]; 
  
  % Need to rotate these into the specified coordinate system  
  % The rotation is not done exactly right, due to not taking into account
  % differences between dbcs and gse. To do in the future
  V = V_dbcs*lmn';
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
  
  MP_grouped_gmm = macroparticle_moments(MP,gm_out.idx);
  pdist_group_gmm = partial_pdist(pdist,MP,gm_out.idx);
  
  % Pressure-based rotation
  % Rotate distribution in the xy-plane so that the pressure is maximized
  % along one dimension. Then 
  [T,vdf] = rotate_vdf_based_on_pressure(pdist,pressure,vel);
  % Fit 1D distributions
  % ...
  
  % Collect output
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


%% Plot results, comparison of kmeans and gmm, only moments
kmN1234_MP = kmN1234_MP_; kmN1234_MP.data = kmN1234_MP.data*0 + 1;
vyN1234_MP = vyN1234_MP_; %vyN1234_MP.data = vyN1234_MP.data*0 + 1;
fig = figure(72);
h = irf_plot(6);

color_gr = [mms_colors('2').^3; mms_colors('2').^0.3; mms_colors('3').^3; mms_colors('3').^0.3];
%color_gr(2) = color_gr(1).^2;
%color_gr(4) = color_gr(3).^2;

lmn_str = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f] (GSE)',lmn');
vxref = -0;
if 0 % n FPI
  hca = irf_panel('n fpi');
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  hca.Title.String = lmn_str;
end
if 1 % n FPI + sum of MP
  hca = irf_panel('n fpi + MP');
  hca.ColorOrder = mms_colors('xyz');
  %c_eval('irf_plot(hca,{ne?,vyN1234_MP_,kmN1234_MP_},''comp'');',ic)
  c_eval('irf_plot(hca,{ne?,vyN1234_MP_*1e6*1.0},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'FPI','vy','k-means'},[0.98 0.98])
  irf_legend(hca,{'FPI','vy'},[0.98 0.98])
  hca.Title.String = lmn_str;
end
if 0 % v FPI = vy MP
  hca = irf_panel('v fpi');
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'v_i','(km/s)'};
end

if 1 % v FPI
  hca = irf_panel('v fpi');
  hca.ColorOrder = mms_colors('xyz1');
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z,vyV_MP.x},''comp'');',ic)
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'v_i','(km/s)'};
end
if 0 % n km MP
  hca = irf_panel('n km MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmN1_MP,kmN2_MP,kmN3_MP,kmN4_MP},'comp');
  hca.YLabel.String = {'n^{km,MP}','(cm^{-3})'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 1 % n vy MP
  hca = irf_panel('n km MP');
  hca.ColorOrder = color_gr;
  irf_plot(hca,{vyN1_MP,vyN2_MP,vyN3_MP,vyN4_MP},'comp');
  hca.YLabel.String = {'n^{vy,MP}','(cm^{-3})'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
end

if 0 % vx km
  hca = irf_panel('vx km');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmV1_MP.x,kmV2_MP.x,kmV3_MP.x,kmV4_MP.x},'comp');
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'v_{ix}^{km,MP}','(km/s)'};
end
if 0 % vy km
  hca = irf_panel('vy km');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmV1_MP.y,kmV2_MP.y,kmV3_MP.y,kmV4_MP.y},'comp');
  hca.YLabel.String = {'v_{iy}^{km,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vz km
  hca = irf_panel('vz km');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmV1_MP.z,kmV2_MP.z,kmV3_MP.z,kmV4_MP.z},'comp');
  hca.YLabel.String = {'v_{iz}^{km,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % n*vx km
  hca = irf_panel('jx km');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmN1_MP*(kmV1_MP.x-vxref)/kmN1234_MP,kmN1_MP*(kmV2_MP.x-vxref)/kmN1234_MP,kmN1_MP*(kmV3_MP.x-vxref)/kmN1234_MP,kmN1_MP*(kmV4_MP.x-vxref)/kmN1234_MP},'comp');
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'nv_{ix}^{km,MP}','(km/s)'};
end
if 0 % n*vy km
  hca = irf_panel('jy km');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmN1_MP*kmV1_MP.y,kmN1_MP*kmV2_MP.y,kmN1_MP*kmV3_MP.y,kmN1_MP*kmV4_MP.y},'comp');
  hca.YLabel.String = {'nv_{iy}^{km,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % n*vz km
  hca = irf_panel('jz km');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{kmN1_MP*kmV1_MP.z,kmN1_MP*kmV2_MP.z,kmN1_MP*kmV3_MP.z,kmN1_MP*kmV4_MP.z},'comp');
  hca.YLabel.String = {'nv_{iz}^{km,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 1 % n*vx vy 1-2
  hca = irf_panel('jx vy');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{vyN1_MP*(vyV1_MP.x-vxref)/vyN1234_MP,vyN2_MP*(vyV2_MP.x-vxref)/vyN1234_MP},'comp');
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'nv_{ix}^{vy,MP}/n_{tot}','(km/s)'};
end
if 1 % n*vy vy 1-2
  hca = irf_panel('jy vy');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{vyN1_MP*vyV1_MP.y/vyN1234_MP,vyN2_MP*vyV2_MP.y/vyN1234_MP},'comp');
  hca.YLabel.String = {'nv_{iy}^{vy,MP}/n_{tot}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
end
if 1 % n*vz vy 1-2
  hca = irf_panel('jz vy');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{vyN1_MP*vyV1_MP.z/vyN1234_MP,vyN2_MP*vyV2_MP.z/vyN1234_MP},'comp');
  hca.YLabel.String = {'nv_{iz}^{vy,MP}/n_{tot}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
end
if 0 % n*vx vy 1-2
  hca = irf_panel('jx vy');
  hca.ColorOrder = color_gr;
  c_eval('toplot? = vyN?_MP*(vyV?_MP.x-vxref)/vyN1234_MP;',1:4)
  irf_plot(hca,{toplot1,toplot2,toplot3,toplot4},'comp');
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'nv_{ix}^{vy,MP}/n_{tot}','(km/s)'};
end
if 0 % n*vy vy 1-2
  hca = irf_panel('jy vy');
  hca.ColorOrder = color_gr;
  c_eval('toplot? = vyN?_MP*(vyV?_MP.y)/vyN1234_MP;',1:4)
  irf_plot(hca,{toplot1,toplot2,toplot3,toplot4},'comp');
  hca.YLabel.String = {'nv_{iy}^{vy,MP}/n_{tot}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
end
if 0 % n*vz vy 1-2
  hca = irf_panel('jz vy');
  hca.ColorOrder = color_gr;
  c_eval('toplot? = vyN?_MP*(vyV?_MP.z)/vyN1234_MP;',1:4)
  irf_plot(hca,{toplot1,toplot2,toplot3,toplot4},'comp');
  hca.YLabel.String = {'nv_{iz}^{vy,MP}/n_{tot}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'v_y<0','vy>0'},[0.98 0.98])
end



if 0 % n gm MP
  hca = irf_panel('n gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmN1_MP,gmN2_MP,gmN3_MP,gmN4_MP},'comp');
  hca.YLabel.String = {'n^{gm,MP}','(cm^{-3})'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vx gm MP
  hca = irf_panel('vx gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1_MP.x,gmV2_MP.x,gmV3_MP.x,gmV4_MP.x},'comp');
  hca.YLabel.String = {'v_{ix}^{gm,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
  hold(hca,'on')
  hca.YLabel.String = 'v_{ix}^{gm,MP} (gm/s)';
  hold(hca,'off')
end
if 0 % vy gm MP
  hca = irf_panel('vy gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1_MP.y,gmV2_MP.y,gmV3_MP.y,gmV4_MP.y},'comp');
  hca.YLabel.String = {'v_{iy}^{gm,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vz gm MP
  hca = irf_panel('vz gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1_MP.z,gmV2_MP.z,gmV3_MP.z,gmV4_MP.z},'comp');
  hca.YLabel.String = {'v_{iz}^{gm,MP}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % n gm model
  hca = irf_panel('n gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmN1,gmN2,gmN3,gmN4},'comp');
  hca.YLabel.String = {'n^{gm,model}','(cm^{-3})'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vx model
  hca = irf_panel('vx gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1.x,gmV2.x,gmV3.x,gmV4.x},'comp');
  hca.YLabel.String = {'v_{ix}^{gm,model}','(km/s)';}
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
  hold(hca,'on')
  hca.YLabel.String = 'v_{ix}^{gm,MP} (gm/s)';
  hold(hca,'off')
end
if 0 % vy model
  hca = irf_panel('vy gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1.y,gmV2.y,gmV3.y,gmV4.y},'comp');
  hca.YLabel.String = {'v_{iy}^{gm,mod}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vz model
  hca = irf_panel('vz gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1.z,gmV2.z,gmV3.z,gmV4.z},'comp');
  hca.YLabel.String = {'v_{iz}^{gm,mod}','(km/s)'};
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % p trace model
  hca = irf_panel('p trace gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmP1.trace/3*1e-6,gmP2.trace/3*1e-6,gmP3.trace/3*1e-6,gmP4.trace/3*1e-6},'comp');
  hca.YLabel.String = 'trace(P)^{gm,mod} (...)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0
  hca = irf_panel('v km');
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{kmV?.x,kmV?.y,kmV?.z},''comp'');',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end
if 0
  hca = irf_panel('n gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmN1,gmN2,gmN3,gmN4},'comp');
  hca.YLabel.String = 'n (cm^{-3})';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0
  hca = irf_panel('n gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmN1_MP,gmN2_MP,gmN3_MP,gmN4_MP},'comp');
  hca.YLabel.String = 'n^{MP} (cm^{-3})';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0
  hca = irf_panel('v fpi');
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
end

irf_zoom(h,'x',times)
irf_zoom(h,'y')
irf_plot_axis_align
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
h(end).XTickLabelRotation = 0;

hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); 
%c_eval('hl(?).LineStyle = ''-'';',1:numel(hl))
c_eval('hl(?).LineWidth = 1.0;',1:numel(hl))

%htmp = findobj(h(7:11),'type','line');
%c_eval('htmp(?).LineStyle = ''none''; htmp(?).Marker = ''*'';',1:numel(htmp))
drawnow
irf_pl_mark(h,time_xline,'black','linestyle','-')
  


%% Plot results, moments and distributions, gaussian mixture model
fig = figure(73);
[h,h2] = initialize_combined_plot('leftright',7,3,6,0.3,'vertical');

lmn_str = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f] (GSE)',lmn');

if 1 % n FPI
  hca = irf_panel('n fpi');
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = 'n (cm^{-3})';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  hca.Title.String = lmn_str;
end
if 1 % v FPI
  hca = irf_panel('v fpi');
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'x','y','z'},[0.98 0.98])
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')
  hold(hca,'off')
  hca.YLabel.String = 'v_i (km/s)';
end
if 1 % n gm MP
  hca = irf_panel('n gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmN1_MP,gmN2_MP,gmN3_MP,gmN4_MP},'comp');
  hca.YLabel.String = 'n^{gm,MP} (cm^{-3})';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 1 % vx gm MP
  hca = irf_panel('vx gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1_MP.x,gmV2_MP.x,gmV3_MP.x,gmV4_MP.x},'comp');
  hca.YLabel.String = 'v_{ix}^{gm,MP} (km/s)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
  hold(hca,'on')   
  irf_plot(hca,irf.ts_scalar(times,-170*ones(nT,1)),'k--')    
  hold(hca,'off')
  hca.YLabel.String = 'v_{ix}^{gm,MP} (gm/s)';
end
if 1 % vy gm MP
  hca = irf_panel('vy gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1_MP.y,gmV2_MP.y,gmV3_MP.y,gmV4_MP.y},'comp');
  hca.YLabel.String = 'v_{iy}^{gm,MP} (km/s)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 1 % vz gm MP
  hca = irf_panel('vz gm MP');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1_MP.z,gmV2_MP.z,gmV3_MP.z,gmV4_MP.z},'comp');
  hca.YLabel.String = 'v_{iz}^{gm,MP} (km/s)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % n gm model
  hca = irf_panel('n gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmN1,gmN2,gmN3,gmN4},'comp');
  hca.YLabel.String = 'n^{gm,model} (cm^{-3})';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vx gm model
  hca = irf_panel('vx gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1.x,gmV2.x,gmV3.x,gmV4.x},'comp');
  hca.YLabel.String = 'v_{ix}^{gm,model} (km/s)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
  hold(hca,'on')
  hca.YLabel.String = 'v_{ix}^{gm,MP} (gm/s)';
  hold(hca,'off')
end
if 0 % vy gm model
  hca = irf_panel('vy gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1.y,gmV2.y,gmV3.y,gmV4.y},'comp');
  hca.YLabel.String = 'v_{iy}^{gm,model} (km/s)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 0 % vz gm model
  hca = irf_panel('vz gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmV1.z,gmV2.z,gmV3.z,gmV4.z},'comp');
  hca.YLabel.String = 'v_{iz}^{gm,model} (km/s)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end
if 1 % p trace gm model
  hca = irf_panel('p trace gm');
  hca.ColorOrder = mms_colors('1234');
  irf_plot(hca,{gmP1.trace/3*1e-6,gmP2.trace/3*1e-6,gmP3.trace/3*1e-6,gmP4.trace/3*1e-6},'comp');
  hca.YLabel.String = 'trace(P)^{gm,model} (...)';
  hca.ColorOrder = mms_colors('1234');
  irf_legend(hca,{'1','2','3','4'},[0.98 0.98])
end

irf_zoom(h,'x',times)
irf_plot_axis_align
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
h(end).XTickLabelRotation = 0;

hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); 
%c_eval('hl(?).LineStyle = ''-'';',1:numel(hl))
c_eval('hl(?).LineWidth = 1.0;',1:numel(hl))

%htmp = findobj(h(7:11),'type','line');
%c_eval('htmp(?).LineStyle = ''none''; htmp(?).Marker = ''*'';',1:numel(htmp))

hmark = irf_pl_mark(h,time_xline,'black','linestyle','-');
  
%%
for it = 1:nT
  pause(0.1)

  hp = findobj(gcf,'type','patch'); delete(hp);
%it = 22;
time = times(it);

time_str = time.utc('MMHHSS_mmm');
print_str = sprintf('ion_dyn_gmm_%s',time_str);


if exist('hmark','var'); delete(hmark); end
c_eval('hmark(?) = irf_pl_mark(h(?),time + nMovMean*0.150*0.5*[-1 1],''k'');',1:numel(h))


group_colors = mms_colors('1234');

% Plot distributions
isub = 1;
if 1 % f(x,y), original total distribution
  hca = h2(isub); isub = isub + 1;  
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  vdf = pdist.reduce('2D',L,M);   
  [h_surf,h_axis,h_all] = vdf.plot_plane(hca);
  h_all.Colorbar.Location = 'south';
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_M (km/s)';
  hold(hca,'on')  
  for iGroup = 1:nGroups
    vx_tmp = out.gm.v{iGroup}(it,1);
    vy_tmp = out.gm.v{iGroup}(it,2);
    vz_tmp = out.gm.v{iGroup}(it,3);
    plot(hca,vx_tmp,vy_tmp,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
  end
  hold(hca,'off')
end
if 1 % f(x,z), original total distribution
  hca = h2(isub); isub = isub + 1;  
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  vdf = pdist.reduce('2D',L,N);   
  vdf.plot_plane(hca);
  [h_surf,h_axis,h_all] = vdf.plot_plane(hca);
  h_all.Colorbar.Location = 'south';
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'v_N (km/s)';
  hold(hca,'on')  
  for iGroup = 1:nGroups
    vx_tmp = out.gm.v{iGroup}(it,1);
    vy_tmp = out.gm.v{iGroup}(it,2);
    vz_tmp = out.gm.v{iGroup}(it,3);
    plot(hca,vx_tmp,vz_tmp,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
  end
  hold(hca,'off')
end
if 1 % f(z,y), original total distribution
  hca = h2(isub); isub = isub + 1;  
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  vdf = pdist.reduce('2D',M,N);     
  [h_surf,h_axis,h_all] = vdf.plot_plane(hca);
  h_all.Colorbar.Location = 'south';
  hca.XLabel.String = 'v_M (km/s)';
  hca.YLabel.String = 'v_N (km/s)';
  hold(hca,'on')  
  for iGroup = 1:nGroups
    vx_tmp = out.gm.v{iGroup}(it,1);
    vy_tmp = out.gm.v{iGroup}(it,2);
    vz_tmp = out.gm.v{iGroup}(it,3);
    plot(hca,vy_tmp,vz_tmp,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);
  end
  hold(hca,'off')
end


% Fitted distributions
vmax = 2400;
nv = 80;
vvec = linspace(-vmax,vmax,nv);
dv = (vvec(2)-vvec(1))*1e3;
[X,Y,Z] = ndgrid(vvec,vvec,vvec);
XYZ = [X(:) Y(:) Z(:)];

for iGroup = 1:nGroups  
  % Make distribtions from gmm
  mu = out.gm.model{it}.mu;
  Sigma = out.gm.model{it}.Sigma;
  Ftmp = out.gm.model{it}.ComponentProportion(iGroup)*mvnpdf(XYZ, mu(iGroup,:), Sigma(:,:,iGroup)); 
  F = reshape(Ftmp,size(X));
  %F = log10(F);

  dolog10 = 0;


  vx_tmp = out.gm.v{iGroup}(it,1);
  vy_tmp = out.gm.v{iGroup}(it,2);
  vz_tmp = out.gm.v{iGroup}(it,3);

  if 1 % f(x,y), original total distribution
    hca = h2(isub); isub = isub + 1;  
    vdf = 1e-2*dv*sum(F,3);
    if dolog10, vdf = log10(vdf); end
    pcolor(hca,vvec,vvec,vdf')    
    shading(hca,'flat')    
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.Box = 'on'; hca.Layer = 'top'; hca.XGrid = 'on'; hca.YGrid = 'on';

    hold(hca,'on')      
    plot(hca,vx_tmp,vy_tmp,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);    
    hold(hca,'off')
  end
  if 1 % f(x,z), original total distribution
    hca = h2(isub); isub = isub + 1;  
    vdf = 1e-2*dv*squeeze(sum(F,2));
    if dolog10, vdf = log10(vdf); end
    pcolor(hca,vvec,vvec,vdf')    
    shading(hca,'flat')    
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    hca.Box = 'on'; hca.Layer = 'top'; hca.XGrid = 'on'; hca.YGrid = 'on';

    hold(hca,'on')      
    plot(hca,vx_tmp,vz_tmp,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);    
    hold(hca,'off')
  end
  if 1 % f(z,y), original total distribution
    hca = h2(isub); isub = isub + 1; 
    vdf = 1e-2*dv*squeeze(sum(F,1));
    if dolog10, vdf = log10(vdf); end
    pcolor(hca,vvec,vvec,vdf')    
    shading(hca,'flat')    
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    hca.Box = 'on'; hca.Layer = 'top'; hca.XGrid = 'on'; hca.YGrid = 'on';

    hold(hca,'on')      
    plot(hca,vy_tmp,vz_tmp,'markeredgecolor','k','markerfacecolor',group_colors(iGroup,:),'Marker','s','linewidth',0.5);    
    hold(hca,'off')
  end
end

if 1 % 1D reduced distributions, observed and gmm, x
  hca = h2(isub); isub = isub + 1;  
  % Distribution from observations
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  vdf_obs = pdist.reduce('1D',L);

  % Make distribtions from gmm
  vdf_model = zeros(nGroups,nv);
  for iGroup = 1:nGroups
    mu = out.gm.model{it}.mu;
    Sigma = out.gm.model{it}.Sigma;
    Ftmp = out.gm.model{it}.ComponentProportion(iGroup)*mvnpdf(XYZ, mu(iGroup,:), Sigma(:,:,iGroup)); 
    F = reshape(Ftmp,size(X));
    vdf_model(iGroup,:) = 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,3),2));
  end        
  
  plot(hca,vdf_obs.depend{1},vdf_obs.data,'color',[0.8 0.6 0.0]);
  %hp = patch(hca,[vdf_obs.depend{1} vdf_obs.depend{1}(1)],[vdf_obs.data vdf_obs.data(1)],1);
  hca.YLabel.String = 'f^{obs} (s/m^3)';   
  hca.ColorOrder = group_colors;
  hold(hca,'on')
  hl = plot(hca,vvec,sum(vdf_model,1)./max(sum(vdf_model,1))*max(vdf_obs.data),'color',[0.6 0.6 0.6],'LineStyle',':');
  hl = plot(hca,vvec,vdf_model./max(sum(vdf_model,1))*max(vdf_obs.data));
  %hl = hl(end:-1:1);
  c_eval('hl(?).Color = group_colors(?,:);',1:nGroups)
  hold(hca,'off')
  hca.XLabel.String = 'v_L (km/s)';
end
if 1 % 1D reduced distributions, observed and gmm, y
  hca = h2(isub); isub = isub + 1;  
  % Distribution from observations
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  vdf_obs = pdist.reduce('1D',M);

  % Make distribtions from gmm
  vdf_model = zeros(nGroups,nv);
  for iGroup = 1:nGroups
    mu = out.gm.model{it}.mu;
    Sigma = out.gm.model{it}.Sigma;
    Ftmp = out.gm.model{it}.ComponentProportion(iGroup)*mvnpdf(XYZ, mu(iGroup,:), Sigma(:,:,iGroup)); 
    F = reshape(Ftmp,size(X));
    vdf_model(iGroup,:) = 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,3),1));
  end        
  
  plot(hca,vdf_obs.depend{1},vdf_obs.data,'color',[0.8 0.6 0.0]);
  %hp = patch(hca,[vdf_obs.depend{1} vdf_obs.depend{1}(1)],[vdf_obs.data vdf_obs.data(1)],1);
  %hp.EdgeColor = 'none';
  hca.YLabel.String = 'f^{obs} (s/m^3)';      
  hca.ColorOrder = group_colors; 
  hold(hca,'on')
  hl = plot(hca,vvec,sum(vdf_model,1)./max(sum(vdf_model,1))*max(vdf_obs.data),'color',[0.6 0.6 0.6],'LineStyle',':');
  hl = plot(hca,vvec,vdf_model./max(sum(vdf_model,1))*max(vdf_obs.data));
  %hl = hl(end:-1:1);
  c_eval('hl(?).Color = group_colors(?,:);',1:nGroups)  
  hold(hca,'off')
  hca.XLabel.String = 'v_M (km/s)';
end
if 1 % 1D reduced distributions, observed and gmm, z
  hca = h2(isub); isub = isub + 1;  
  % Distribution from observations
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  vdf_obs = pdist.reduce('1D',N);

  % Make distribtions from gmm
  vdf_model = zeros(nGroups,nv);
  for iGroup = 1:nGroups
    mu = out.gm.model{it}.mu;
    Sigma = out.gm.model{it}.Sigma;
    Ftmp = out.gm.model{it}.ComponentProportion(iGroup)*mvnpdf(XYZ, mu(iGroup,:), Sigma(:,:,iGroup)); 
    F = reshape(Ftmp,size(X));
    vdf_model(iGroup,:) = 1e-2*1e-2*dv*dv*squeeze(sum(sum(F,2),1));
  end        
  
  plot(hca,vdf_obs.depend{1},vdf_obs.data,'color',[0.8 0.6 0.0]);
  %hp = patch(hca,[vdf_obs.depend{1} vdf_obs.depend{1}(1)],[vdf_obs.data vdf_obs.data(1)],1);
  hca.YLabel.String = 'f^{obs} (s/m^3)';   
  hca.ColorOrder = group_colors;
  hold(hca,'on')
  hl = plot(hca,vvec,sum(vdf_model,1)./max(sum(vdf_model,1))*max(vdf_obs.data),'color',[0.6 0.6 0.6],'LineStyle',':');
  hl = plot(hca,vvec,vdf_model./max(sum(vdf_model,1))*max(vdf_obs.data));
  %hl = hl(end:-1:1);
  c_eval('hl(?).Color = group_colors(?,:);',1:nGroups)   
  hold(hca,'off')
  hca.XLabel.String = 'v_N (km/s)';
end


c_eval('axis(h2(?),''square'')',1:numel(h2))
%hlinks = linkprop(h2,{'XLim','YLim','CLim'});
hlinks = linkprop(h2(1:15),{'XLim','YLim'});
hlinks_obs = linkprop(h2(1:3),{'CLim'});
hlinks_model = linkprop(h2(4:15),{'CLim'});
%colormap(flipdim(pic_colors('thermal'),1))
colormap(pic_colors('candy_gray'))


hb = findobj(gcf,'type','colorbar');
c_eval('hb(?).Position(1) = hb(?).Position(1) + 0.6*hb(?).Position(3);',1:numel(hb))
c_eval('hb(?).Position(3) = 0.4*hb(?).Position(3);',1:numel(hb))
c_eval('hb(?).Position(4) = 0.01;',1:numel(hb))
hl = findobj(gcf,'type','line'); hl = hl(end:-1:1); 
c_eval('hl(?).LineWidth = 1 ;',1:numel(hl))

cn.print(print_str,'path',[printpath 'gmm_non-diag_cov' filesep])

end % end time loop


















