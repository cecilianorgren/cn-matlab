nMovMean = 5; % Number of distributions for the moving average
elim = [200 Inf]; % Energy limit
N_MP = 100000; % Number of macroparticles
nGroupsMax = 4;  % number of classes/groups for kmeans and gmm

T = nMovMean*0.150;
pdist_all = PD_clean.movmean(nMovMean);
it = 1;pdist_all.length-10;
pdist = pdist_all(it);
time = pdist.time;

c_eval('scpot = mean(scPot?.tlim(time + T* 0.15*0.5*[-1 1]).data,1);',ic)
scpot = irf.ts_scalar(time,scpot);
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


aic = zeros(nGroupsMax,1);
bic = zeros(nGroupsMax,1);
nLogL = zeros(nGroupsMax,1);
for nGroups = 1:nGroupsMax  
  initial_guess = [0 -1000 -1000;...
                   0 -1000 +1000;...
                   0 +1000 -1000;...
                   0 +1000 +1000;...
                   0  0000  0000;...
                   0  0000  0000];
  initial_guess = initial_guess(1:nGroups,:);
  S.mu = initial_guess; % same as for kmeans         
  S.Sigma = repmat([500 10 10; 10 500 10; 10 10 500],[1 1 nGroups]);  
  S.ComponentProportion = repmat(1,[1,nGroups]);

  X = [MP.vx, MP.vy, MP.vz];
  gm = fitgmdist(X,nGroups,'Start',S,'SharedCovariance',false);
  aic(nGroups) = gm.AIC;
  bic(nGroups) = gm.BIC;
  nLogL(nGroups) = gm.NegativeLogLikelihood;
end

h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,1:nGroupsMax,aic/nMP)
hca.XLabel.String = 'K';
hca.YLabel.String = 'AIC';

hca = h(isub); isub = isub + 1;
plot(hca,1:nGroupsMax,bic/nMP)
hca.XLabel.String = 'K';
hca.YLabel.String = 'BIC';

hca = h(isub); isub = isub + 1;
plot(hca,1:nGroupsMax,nLogL)
hca.XLabel.String = 'K';
hca.YLabel.String = 'nlogL';
