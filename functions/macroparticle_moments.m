function varargout = macroparticle_moments(MP,idx)
% MACROPARTICLE_MOMENTS Calculates moments based on group index.

fieldnames = fields(MP);
nGroups = numel(unique(idx));

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

varargout{1} = MP_grouped;