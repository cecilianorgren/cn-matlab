function varargout = gmm_evaluate_f(gm,XYZ,ntot,varargin)

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

%vx = vxyz(:,1);
%[X,Y,Z] = ndgrid(vx,vy,vz);
%XYZ = [X(:) Y(:) Z(:)];
%dvx = vx(2)-vx(1);
%dvy = vy(2)-vy(1);
%dvz = vz(2)-vz(1);

mu = gm.mu; % km/s
Sigma = gm.Sigma; % (km/s)^2
compProp = gm.ComponentProportion;
numComp = gm.NumComponents;
[np,ndim] = size(XYZ);
Ftot = zeros([np 1]);
Fcomp = zeros([np numComp]);

for iComp = 1:gm.NumComponents             
  ntot_km3 = ntot;%ntot*1e15; % cm^-3 -> km^-3
  Ftmp = ntot_km3*compProp(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); % s^3/km^6      
%  Ftmp = reshape(Ftmp,size(X));  
  Ftot = Ftot + Ftmp;
  Fcomp(:,iComp)  = Ftmp;
  % sum(Ftot(:).*Fobs.dv(:)) = km^-3
end

if doGroup
  nGroups = numel(groups);
  Fgrouped = zeros([np nGroups]);
  for iGroup = 1:nGroups
    group_tmp = groups{iGroup};
    Fgrouped(:,iGroup) = sum(Fcomp(:,group_tmp),2);    
  end
end


switch nargout
  case 1
    varargout{1} = Ftot;
  case 2
    varargout{1} = Ftot;
    varargout{2} = Fcomp;
  case 3
    varargout{1} = Ftot;
    varargout{2} = Fcomp;
    varargout{3} = Fgrouped;
end