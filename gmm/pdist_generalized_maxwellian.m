function PD = pdist_generalized_maxwellian(pd,w,mu,sigma)
% pdist_generalized_maxwellian from weights, mu, and Sigma
%   PD = pdist_generalized_maxwellian(grid,w,mu,Sigma)
%
%     grid - grid to evaluate f on, options:
%            - PDist, [vx,vy,vz] = PDist.v; - output is a PDist
%            - (NOT IMPLEMENTED) numeric, R = [Nx3], where vx = R(:,1); vy = R(:,2); vz = R(:,3); 
%     w - [K,nt], weight or density
%     mu - [K,3,nt], mean/velocity
%     Sigma - [3,3,K,nt], covariance matrix/v_th^2 (strange ordering is due to gmdistribution output format)

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
    
    % Velocity in rest frame
    dvx = vx(it,:,:,:)-mu(iK,1,it);
    dvy = vy(it,:,:,:)-mu(iK,2,it);
    dvz = vz(it,:,:,:)-mu(iK,3,it);
  
    dV = [dvx(:), dvy(:), dvz(:)];
    
    % Exponent
    y = dV / L;
    q = sum(y.^2, 2); % L2 norm
    
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