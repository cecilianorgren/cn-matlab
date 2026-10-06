function PD = generalized_maxwellian_macroparticle(pd,w,mu,sigma,N)

units = irf_units;

% Convert to SI units
pd = pd.convertto('s^3/m^6');
[vx,vy,vz] = pd.v; % velocity coordinates at center energy of bin 

nt = length(pd.time);
if isa(w,'TSeries')
  nK = 1;
end
if isa(w,'TSeries')
  w = w.data;
  w = reshape(w,[nK,nt]);
  w = w*1e-3; % what's this unit, should not be cm^-3 -> km^-3
end
if isa(mu,'TSeries')
  mu_inp = permute(mu.data,[2 1]);
  mu = reshape(mu_inp,[nK,3,nt]);
end
if isa(sigma,'TSeries') % input is T
  T = permute(sigma.data,[2 3 1 ]); % The strange ordering is because the gmdistribution returns the sigma like this
  T = reshape(T,[size(T,[1 2]) 1 size(T,3)]);
  sigma_m2s2 = 2*T*units.eV/units.mp; % v^2 [(m/s)^2] ... mv^2/2=kBT
  sigma_km2s2 = sigma_m2s2/1e6;
  sigma = sigma_km2s2;
end

nK = size(w,1);

% Construct generalized Maxwellian
counts_all = zeros(size(vx));
for it = 1:nt
  for iK = 1:nK
    Nk = N*w(iK,it)/sum(w(:,it)); % weight/density

    % One can go this way instead of taking the determinant of the covariance matrix
    % Sigma = LL^T
    % (v-mu)^T * Sigma^-1 * (v-mu) = ||(v-mu)*L^-T||^2
    L = chol(sigma(:,:,iK,it), 'lower'); 

    v = mu(iK,:,it) + randn(Nk,3)*L.';
    v = -v; % where they're coming from, not where they're going

    % Bin to FPI grid
    vabs = sqrt(sum(v.^2,2));
    E = pd.mass*vabs.^2*1e6/(2*units.eV);
    azim = atan2d(v(:,2),v(:,1))+180;
    pol = acosd(v(:,3)./vabs);

    grid.E = [pd.depend{1}(it,:)-pd.ancillary.delta_energy_minus(it,:) pd.depend{1}(it,end)+pd.ancillary.delta_energy_plus(it,end)];
    daz = diff(pd.depend{2},1,2);
    grid.az = [pd.depend{2}(it,1:end-1)-0.5*daz(it,:) pd.depend{2}(it,end)-0.5*daz(it,end) pd.depend{2}(it,end)+0.5*daz(it,end)]; 
    dpol = diff(pd.depend{3},1,2);
    grid.pol = [pd.depend{3}(:,1:end-1)-0.5*dpol pd.depend{3}(:,end)-0.5*dpol(:,end) pd.depend{3}(:,end)+0.5*dpol(:,end)]; 

    [counts] = histcn([E azim,pol],grid.E,grid.az,grid.pol);

    % Gather results
    counts = reshape(counts, [1 size(vx,[2 3 4])]);
    counts_all(it,:,:,:) = counts_all(it,:,:,:) + counts;
  end
end
f = counts_all./pd.d3v('mat')/N;

% Make PDist file for output
PD = pd;
PD.data = f;
PD = PD.convertto('s^3/km^6');
PD.name = 'generalized_maxwellian';