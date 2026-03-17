function varargout = gmm_moments_from_macroparticles(MP,idx,igr)
% MACROPARTICLE_MOMENTS Calculates moments based on group index.
units = irf_units;
fieldnames = fields(MP);
nGroups = numel(unique(igr));

idx_tmp = [];
for iGroup = 1:nGroups
  idx_tmp = cat(1,idx_tmp,find(idx==(igr(iGroup))));
end
nMP = numel(MP.iDep1);
%for iGroup = 1:nGroups
  %idx_tmp = find(idx==iGroup);   
  %idx_tmp = idx;
  for iField = 1:numel(fieldnames)
    var = MP.(fieldnames{iField});
    if numel(var) == nMP
      MPtmp.(fieldnames{iField}) = var(idx_tmp);
    end
  end
  
  % Calculate speed and density
  n = sum(MPtmp.dn);
  jx = sum(MPtmp.vx.*MPtmp.dn);
  jy = sum(MPtmp.vy.*MPtmp.dn);
  jz = sum(MPtmp.vz.*MPtmp.dn);
  vx = jx/n;
  vy = jy/n;
  vz = jz/n;
  dpxx = MPtmp.dn.*(MPtmp.vx-vx).*(MPtmp.vx-vx);
  dpyy = MPtmp.dn.*(MPtmp.vy-vy).*(MPtmp.vy-vy);
  dpzz = MPtmp.dn.*(MPtmp.vz-vz).*(MPtmp.vz-vz);
  dpxy = MPtmp.dn.*(MPtmp.vx-vx).*(MPtmp.vy-vy);
  dpxz = MPtmp.dn.*(MPtmp.vx-vx).*(MPtmp.vz-vz);
  dpyz = MPtmp.dn.*(MPtmp.vy-vy).*(MPtmp.vz-vz);
  pxx = sum(dpxx);
  pyy = sum(dpyy);
  pzz = sum(dpzz);
  pxy = sum(dpxy);
  pxz = sum(dpxz);
  pyz = sum(dpyz);


  MPtmp.n = n;
  MPtmp.jx = jx;
  MPtmp.jy = jy;
  MPtmp.jz = jz;
  MPtmp.vx = vx;
  MPtmp.vy = vy;
  MPtmp.vz = vz;
  MPtmp.pxx = pxx;
  MPtmp.pyy = pyy;
  MPtmp.pzz = pzz;
  MPtmp.pxy = pxy;
  MPtmp.pxz = pxz;
  MPtmp.pyz = pyz;

  %MP_grouped{iGroup} = MPtmp;
  MP_grouped = MPtmp;
%end

varargout{1} = MP_grouped;