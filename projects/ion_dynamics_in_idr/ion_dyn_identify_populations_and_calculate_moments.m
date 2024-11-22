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

c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',1:4);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)

% Local coordinate system
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

lmn = lmn_gse;
lmn = lmn_edr;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);

c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)

%% Run through distributions
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');

nMovMean = 5; % Number of distributions for the moving average
elim = [200 Inf]; % Energy limit
N = 100000; % Number of macroparticles
nGroups = 4;  % number of classes/groups for kmeans and gmm

% Prepare the distribution
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)

for dt = -5
  time = time_xline + dt;  
  %pdist = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).tlim(time + 0.15*0.5*[-1 1]).elim([300 Inf]);
  pdist = pdist_all.tlim(time + 0.15*0.5*[-1 1]);
  
  % These should be in dbcs
  vel = mean(gseVi3.tlim(time + nMovMean*0.15*0.5*[-1 1]).data,1);
  pressure = mean(gsePi3.tlim(time +nMovMean* 0.15*0.5*[-1 1]).data,1);

  % Initial macroparticles
  MP = pdist.elim([50 Inf]).macroparticles('ntot',N,'skipzero',1);
  V = [MP.vx, MP.vy, MP.vz]; 

  % Initial guess for 4 populations
  initial_guess = [0 -1000 -1000;...
                   0 -1000 +1000;...
                   0 +1000 -1000;...
                   0 +1000 +1000];
  % kmeans  
  [idx, c, sumd, d] = kmeans(V, nGroups,'Start',initial_guess);
  %s = silhouette(V, idx, 'sqeuclid');
  
  
  MP_grouped_km = macroparticle_moments(MP,idx);
  pdist_group_km = partial_pdist(pdist,MP,idx);
  
  % Gaussian mixture model
  S.mu = initial_guess;
  
  %S.Sigma = [500 500 500]'; % this is not right format
  S.Sigma = 500*ones(3,3,nGroups);
  S.Sigma = 500*ones(1,3,nGroups);
  S.ComponentProportion = [1 1 1 1];
  
  %gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],4,'Start',S);
  gm_out = gaussian_mixture_model([MP.vx, MP.vy, MP.vz],4,'Start',S,'CovType',...
    'diagonal');
  
  MP_grouped_gmm = macroparticle_moments(MP,gm_out.idx);
  pdist_group_gmm = partial_pdist(pdist,MP,gm_out.idx);
  
  % Pressure-based rotation
  % Rotate distribution in the xy-plane so that the pressure is maximized
  % along one dimension. Then 
  [T,vdf] = rotate_vdf_based_on_pressure(pdist,pressure,vel);
  % Fit 1D distributions
  % ...
  
  % Plot comparison
  fig = figure(71);
  h = setup_subplots(1,1);
  
end















