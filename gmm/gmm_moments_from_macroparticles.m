function varargout = gmm_moments_from_macroparticles(MP,gm,igr,varargin)
% MACROPARTICLE_MOMENTS Calculates moments based on group index.
%
% By default, cluster labels are remapped to match the component sorting used
% elsewhere (see `gmm_sort`): sorted component index 1 corresponds to the
% smallest temperature proxy. This makes `igr` stable in time when you want
% to compare macroparticle moments with GMM moments computed using `isort`.
%
% Options (name-value):
%   'sort' (default true): alias for DoSortComponents, if true remap cluster
%     labels to sorted component indices (coldest=1, ...).
%   'DoSortComponents' (default true): if true, remap cluster labels to
%     sorted component indices (coldest=1, ...).
doSortComponents = true;
if ~isempty(varargin)
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'sort'
        doSortComponents = varargin{iarg+1};
      case 'dosortcomponents'
        doSortComponents = varargin{iarg+1};
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end

units = irf_units;
fieldnames = fields(MP);
nGroups = numel(igr);
nt = numel(MP);
nGM = size(gm,2);
moms = cell(nt,nGM,nGroups);
isort_used = cell(nt,nGM);

for it = 1:nt
  MPtmp = MP(it);  
  nMP = numel(MPtmp.iDep1);
  R = [MPtmp.vx, MPtmp.vy, MPtmp.vz];
  for iK = 1:nGM
    idx = gm{it,iK}.cluster(R);
    if doSortComponents
      isort = gmm_sort(gm{it,iK});
      inv_isort = zeros(1,numel(isort));
      inv_isort(isort) = 1:numel(isort);
      idx = inv_isort(idx);
      isort_used{it,iK} = isort;
    else
      isort_used{it,iK} = 1:gm{it,iK}.NumComponents;
    end
    for iGroup = 1:nGroups % Get the indices of all the members of the group
      MPgroup = MPtmp;
      idx_group = [];
      group = igr{iGroup};
      for iSubGroup = 1:numel(group)
        idx_group = cat(1,idx_group,find(idx==group(iSubGroup)));
      end    
     % idx_group(1:10)
    
      for iField = 1:numel(fieldnames) % Remove all indices that are not belonging to the group
        var = MPtmp.(fieldnames{iField});
        if isnumeric(var) && numel(var) == nMP
          MPgroup.(fieldnames{iField}) = var(idx_group);
        end
      end
      %MPgroup

      % Calculate speed and density
      MPgroup.dn  = MPgroup.df.*MPgroup.dv;
      n = sum(MPgroup.dn);
      jx = sum(MPgroup.vx.*MPgroup.dn);
      jy = sum(MPgroup.vy.*MPgroup.dn);
      jz = sum(MPgroup.vz.*MPgroup.dn);
      vx = jx/n;
      vy = jy/n;
      vz = jz/n;
      dpxx = MPgroup.dn.*(MPgroup.vx-vx).*(MPgroup.vx-vx);
      dpyy = MPgroup.dn.*(MPgroup.vy-vy).*(MPgroup.vy-vy);
      dpzz = MPgroup.dn.*(MPgroup.vz-vz).*(MPgroup.vz-vz);
      dpxy = MPgroup.dn.*(MPgroup.vx-vx).*(MPgroup.vy-vy);
      dpxz = MPgroup.dn.*(MPgroup.vx-vx).*(MPgroup.vz-vz);
      dpyz = MPgroup.dn.*(MPgroup.vy-vy).*(MPgroup.vz-vz);
      pxx = sum(dpxx);
      pyy = sum(dpyy);
      pzz = sum(dpzz);
      pxy = sum(dpxy);
      pxz = sum(dpxz);
      pyz = sum(dpyz);
    
      moms{it,iK,iGroup}.n = n;
      moms{it,iK,iGroup}.jx = jx;
      moms{it,iK,iGroup}.jy = jy;
      moms{it,iK,iGroup}.jz = jz;
      moms{it,iK,iGroup}.vx = vx;
      moms{it,iK,iGroup}.vy = vy;
      moms{it,iK,iGroup}.vz = vz;
      moms{it,iK,iGroup}.pxx = pxx;
      moms{it,iK,iGroup}.pyy = pyy;
      moms{it,iK,iGroup}.pzz = pzz;
      moms{it,iK,iGroup}.pxy = pxy;
      moms{it,iK,iGroup}.pxz = pxz;
      moms{it,iK,iGroup}.pyz = pyz;
      % MPtmp.n = n;
      % MPtmp.jx = jx;
      % MPtmp.jy = jy;
      % MPtmp.jz = jz;
      % MPtmp.vx = vx;
      % MPtmp.vy = vy;
      % MPtmp.vz = vz;
      % MPtmp.pxx = pxx;
      % MPtmp.pyy = pyy;
      % MPtmp.pzz = pzz;
      % MPtmp.pxy = pxy;
      % MPtmp.pxz = pxz;
      % MPtmp.pyz = pyz;
    end
  end
end

%nMP = numel(MP.iDep1);
%for iGroup = 1:nGroups
  %idx_tmp = find(idx==iGroup);   
  %idx_tmp = idx;

  
  

  %MP_grouped{iGroup} = MPtmp;
  %MP_grouped = MPtmp;
%end

varargout{1} = moms;
if nargout > 1
  varargout{2} = isort_used;
end