function PD = generalized_maxwellian_macroparticle(pd,w,mu,sigma,N)

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
  w = w*1e-3;
end
if isa(mu,'TSeries')
  mu_inp = permute(mu.data,[2 1]);
  mu = reshape(mu_inp,[nK,3,nt]);
end
if isa(sigma,'TSeries') % input is T
  T = permute(sigma.data,[2 3 1 ]); % The strange ordering is because the gmdistribution returns the sigma like this
  T = reshape(T,[size(T,[1 2]) 1 size(T,3)]);
  units = irf_units;
  sigma_m2s2 = T*units.eV/units.mp; % v^2/2 [(m/s)^2]
  sigma_km2s2 = sigma_m2s2/1e6;
  sigma = sigma_km2s2;
end

nK = size(w,1);

% Construct generalized Maxwellian
f_all = zeros(size(vx));
for it = 1:nt
  for iK = 1:nK
    n = w(iK,it); % weight/density

    % One can go this way instead of taking the determinant of the covariance matrix
    % Sigma = LL^T
    % (v-mu)^T * Sigma^-1 * (v-mu) = ||(v-mu)*L^-T||^2
    L = chol(sigma(:,:,iK,it), 'lower'); 

    v = mu(iK,1,it) + randn(N,3)*L.';
    r = sqrt(sum(v.^2,2));
    azim = atan2d(v(:,2),v(:,1));
    pol = acosd(v(:,3)./r);

    ...
    ...
    ...
    ...
    ...

    % Maxwellian
    f = n / ((2*pi)^(3/2) * prod(diag(L))) .* exp(-0.5*q); 
    f = reshape(f, size(dvx));

    % Gather results
    f_all(it,:,:,:) = f_all(it,:,:,:) + f;
  end
end

% Make PDist file for output
PD = pd;
PD.data = f_all;
PD = PD.convertto('s^3/km^6');
PD.name = 'generalized_maxwellian';