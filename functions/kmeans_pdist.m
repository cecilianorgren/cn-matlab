function varargout = kmeans_pdist(pdist,nGroups,nMP,varargin)
% KMEANS_PDIST(PDist,nGroups,nMacoParticles)

% Make macroparticles from pdist
MP = pdist.macroparticles('ntot',nMP,'skipzero',1);
V = [MP.vx, MP.vy, MP.vz];

% Do kmeans clustering
[idx, c, sumd, d] = kmeans(V, nGroups, varargin{:});

% Group macro particles by k-group
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

% Make nGroups PDists based on the dominant class for each bin
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

out.t = pdist.time;
out.MP = MP;
out.MP_grouped = MP_grouped;
out.pdist = pdist;
out.pdist_grouped = pdist_group;

varargout{1} = out;
