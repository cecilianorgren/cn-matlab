function varargout = gmm_get_F_timeseries(gm,vx,vy,vz,ntot,varargin)

isortOpt = [];
doGroup = 0;
if ~isempty(varargin)
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'isort'
        isortOpt = varargin{iarg+1};
      case 'group'
        groups = varargin{iarg+1};
        doGroup = 1;
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end


[X,Y,Z] = ndgrid(vx,vy,vz);
XYZ = [X(:) Y(:) Z(:)];
%dvx = vx(2)-vx(1);
%dvy = vy(2)-vy(1);
%dvz = vz(2)-vz(1);

if isa(gm,'gmdistribution')
  gm_cell = {gm};
else isa(gm,'cell')
  gm_cell = gm;
end

if doGroup && isa('groups','numeric')
  groups = {groups};
end

[nt,nK] = size(gm_cell);
Fcomp_cell = cell(nt,nK);
Fgrouped_cell = cell(nt,nK);

for it = 1:nt
  gm = gm_cell;
  mu = gm.mu; % km/s
  Sigma = gm.Sigma; % (km/s)^2
  compProp = gm.ComponentProportion;
  numComp = gm.NumComponents;
  Ftot = zeros(size(X));
  Fcomp = zeros([size(X) numComp]);
  
  for iComp = 1:gm.NumComponents             
    ntot_km3 = ntot*1e15; % cm^-3 -> km^-3
    Ftmp = ntot_km3*compProp(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); % s^3/km^6      
    Ftmp = reshape(Ftmp,size(X));  
    Ftot = Ftot + Ftmp;
    Fcomp(:,:,:,iComp)  = Ftmp;
    % sum(Ftot(:).*Fobs.dv(:)) = km^-3
  end
  
  if doGroup
    nGroups = numel(groups{it,iK});
    Fgrouped = zeros([size(X) nGroups]);
    for iGroup = 1:nGroups
      group_tmp = groups{it,iK,iGroup}; % this is not ok probably
      Fgrouped(:,:,:,iGroup) = sum(Fcomp(:,:,:,group_tmp),4);    
    end
  end
  Ftot_cell{it,iK} = Ftot;
  Fcomp_cell{it,iK} = Fcomp;
  Fgrouped_cell{it,iK} = Fgrouped;
end

switch nargout
  case 1
    varargout{1} = Ftot_cell;
  case 2
    varargout{1} = Ftot_cell;
    varargout{2} = Fcomp_cell;
  case 3
    varargout{1} = Ftot_cell;
    varargout{2} = Fcomp_cell;
    varargout{3} = Fgrouped_cell;
end