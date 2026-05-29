function varargout = total_covariance(w0,mu0,Sigma0,ksub)

% Sub components to include
if not(exist('kvec','var'))
  ksub = 1:numel(w0); % if not give, merge all
end

% Number of components
K = numel(ksub); 

% Total weight
wS = sum(w0(ksub));

% Normalized weight
wSnorm = w0(ksub)/wS;

% Normalized mean
muS = (w0(ksub) * mu0(ksub,:)) / wS; % 1x3

% Total covariance (within + between)
SigS = zeros(3,3); % Initialize covariance matrix

for ii = 1:K % sum over components
  k = ksub(ii);
  dmu = (mu0(k,:) - muS)';
  SigS = SigS + w0(k) * (Sigma0(:,:,k) + (dmu*dmu')); % law of shared covariances
end
SigS = SigS / wS;

if nargout == 1
  out.w = wS;
  out.mu = muS;
  out.Sigma = SigS;
  varargout{1} = out;
elseif nargout == 3
  varargout{1} = wS;
  varargout{2} = muS;
  varargout{3} = SigS;
end
  

