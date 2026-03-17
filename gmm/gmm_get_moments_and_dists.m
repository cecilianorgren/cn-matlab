function varargout = gmm_get_moments_and_dists(times,gm,xmid,ymid,zmid,ntot,varargin)

units = irf_units;

% Extract info
nt = size(gm,1);
K = gm{1}.NumComponents;

% Backward compatibility + options
% Old signature: (..., ntot, isort)
% New usage: (..., ntot, 'sort', true/false)
doSort = false;
isort_fixed = [];
if ~isempty(varargin)
  if isnumeric(varargin{1})
    isort_fixed = varargin{1};
    varargin = varargin(2:end);
  end
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'sort'
        doSort = varargin{iarg+1};
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end
if isempty(isort_fixed)
  isort_fixed = 1:K;
end

% Set up grid for vdf computation
dvx = xmid(2)-xmid(1); 
dvy = ymid(2)-ymid(1);
dvz = zmid(2)-zmid(1);
nvx = numel(xmid);
nvy = numel(ymid);
nvz = numel(zmid);
[X,Y,Z] = ndgrid(xmid,ymid,zmid);
XYZ = [X(:) Y(:) Z(:)];

% Initalize variables
T_tens = zeros(nt,3,3,K);
vx = zeros(nt,K);
vy = zeros(nt,K);
vz = zeros(nt,K);
fx = zeros(nt,nvx,K);
fy = zeros(nt,nvy,K);
fz = zeros(nt,nvz,K);
n = zeros(nt,K);
Ftot = zeros(size(X));

for ii = 1:nt
  if doSort
    isort = gmm_sort(gm{ii});
  else
    isort = isort_fixed;
  end
  mu = gm{ii}.mu(isort,:);
  Sigma = gm{ii}.Sigma(:,:,isort);
  compProp = gm{ii}.ComponentProportion(isort);
  for iComp = 1:K     
    % Reduced distributions
    Ftmp = compProp(iComp)*mvnpdf(XYZ, mu(iComp,:), Sigma(:,:,iComp)); 
    Ftmp = Ftmp;
    F{iComp} = reshape(Ftmp,size(X));  
    Ftot = Ftot + F{iComp};
    fx(ii,:,iComp) = sum(F{iComp},[2 3])*ntot(ii)*dvy*dvz*1e3;
    fy(ii,:,iComp) = sum(F{iComp},[1 3])*ntot(ii)*dvx*dvz*1e3;
    fz(ii,:,iComp) = sum(F{iComp},[1 2])*ntot(ii)*dvx*dvy*1e3;

    % Moments
    % Velocity
    vx(ii,iComp) = mu(iComp,1);
    vy(ii,iComp) = mu(iComp,2);
    vz(ii,iComp) = mu(iComp,3);
    % Density
    n(ii,iComp) = compProp(iComp)*ntot(ii);
    % Temperature tensor
    T_tens(ii,:,:,iComp) = Sigma(:,:,iComp)*1e6*units.mp/2/units.eV;
  end
end

% Sum over components
fx_tot =  sum(fx,3);
fy_tot =  sum(fy,3);
fz_tot =  sum(fz,3);

% Construct TSeries
% Component-wise
for iComp = 1:K
  tsFx{iComp} = PDist(times,fx(:,:,iComp),'1Dcart',xmid); tsFx{iComp}.units = 's/m^4';
  tsFy{iComp} = PDist(times,fy(:,:,iComp),'1Dcart',ymid); tsFy{iComp}.units = 's/m^4';
  tsFz{iComp} = PDist(times,fz(:,:,iComp),'1Dcart',zmid); tsFz{iComp}.units = 's/m^4';
  tsT{iComp} = irf.ts_tensor_xyz(times, T_tens(:,:,:,iComp));
  tsV{iComp} = irf.ts_vec_xyz(times, [vx(:,iComp),vy(:,iComp),vz(:,iComp)]);
  tsN{iComp} = irf.ts_scalar(times, n(:,iComp));
end

% Summed together, maybe this is not necessary?
tsFx_tot = PDist(times,fx_tot,'1Dcart',xmid);
tsFy_tot = PDist(times,fy_tot,'1Dcart',ymid);
tsFz_tot = PDist(times,fz_tot,'1Dcart',zmid);
n_vec = cellfun(@(x) x.data, tsN, 'UniformOutput', false);
tsN_tot = irf.ts_scalar(times,sum(cat(2,n_vec{:}),2));

if nargout == 3
  varargout = {tsN,tsV,tsT};
elseif nargout == 6
  varargout = {tsN,tsV,tsT,tsFx,tsFy,tsFz};
elseif nargout == 10
  varargout = {tsN,tsV,tsT,tsFx,tsFy,tsFz,tsN_tot,tsFx_tot,tsFy_tot,tsFz_tot};
end