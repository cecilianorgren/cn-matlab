% Compare the temprature of cold population with the instrument bin width,
% to see if any apparent heating is only due to the increasing bin widths
% at higher speeds
R = [MP.vx, MP.vy, MP.vz];

K = 2;
for it = 1:nt
  time = times(it);
  pdist = PD.tlim(time+0.5*0.150*[-1 1]);
  %nMP = nansum(round(pdist.data(~isnan(pdist.data)))); % total number of counts
  
  scpot = mean(scPot.tlim(time + 0.5*0.15*[-1 1]).data,1);
  scpot = irf.ts_scalar(time,scpot);
     
   MP = pdist.macroparticles('ntot',nMP,'skipzero',1,'scpot',scpot);
   nMP = numel(MP.dv);
   MP.dn = MP.df.*MP.dv;
   V_dbcs = [MP.vx, MP.vy, MP.vz]; 
  
   % Need to rotate these into the specified coordinate system  
   % The rotation is not done exactly right, due to not taking into account
   % differences between dbcs and gse. To do in the future
   %V = V_dbcs*lmn';
   V = V_dbcs;
   MP.vx = V(:,1);
   MP.vy = V(:,2);
   MP.vz = V(:,3);

   idx = cluster(gm{it,K},R);
   ebin = MP.iDep1(idx==1);
   dEbin = iPDist.ancillary.delta_energy_plus - iPDist.ancillary.delta_energy_minus;
   T = tsT{1}.trace/3;