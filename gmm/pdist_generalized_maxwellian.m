function PD = pdist_generalized_maxwellian(pd,w,mu,sigma)
% pdist_generalized_maxwellian from weights, mu, and Sigma
%   PD = pdist_generalized_maxwellian(pd,w,mu,sigma)
%
%     pd - PDist, used for velocity grid and output format
%     w - [nt,1], weight or density
%     mu - [nt,3], mean/velocity
%     sigma - [3,3,nt], covariance matrix/v_th^2

% Convert to SI units
pd = pd.convertto('s^3/m^6');
[vx,vy,vz] = pd.v;

% Construct generalized Maxwellian
f_all = zeros(size(vx));
for ii = 1:length(pd.time)
  n = w(ii); % weight/density

  % One can go this way instead of taking the determinant of the covariance matrix
  % Sigma = LL^T
  % (v-mu)^T * Sigma^-1 * (v-mu) = ||(v-mu)*L^-T||^2
  L = chol(sigma(:,:,ii), 'lower'); 
  
  % Velocity in rest frame
  dvx = vx(ii,:,:,:)-mu(ii,1);
  dvy = vy(ii,:,:,:)-mu(ii,2);
  dvz = vz(ii,:,:,:)-mu(ii,3);

  dV = [dvx(:), dvy(:), dvz(:)];
  
  % Exponent
  y = dV / L;
  q = sum(y.^2, 2); % L2 norm
  
  % Maxwellian
  f = n / ((2*pi)^(3/2) * prod(diag(L))) .* exp(-0.5*q); 
  f = reshape(f, size(dvx));

  % Gather results
  f_all(ii,:,:,:) = f;
end

% Make PDist file for output
PD = pd;
PD.data = f_all;
PD = PD.convertto('s^3/km^6');