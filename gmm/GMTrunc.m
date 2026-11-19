classdef GMTrunc
  % GMTRUNC Gaussian mixture observed through a speed shell, vmin < |v| < vmax.
  %
  %   The macroparticles only exist inside the energy range of the instrument,
  %   a spherical shell of velocity space (in the frame of the instrument). A
  %   fast beam can be cut off by the upper limit. The components of a GMTrunc are
  %   the Gaussians of the UNtruncated mixture: the observed density is
  %
  %     p_obs(v) = sum_k pi_k N(v; mu_k, Sigma_k) / A   inside the shell, 0 outside,
  %     A = sum_k pi_k Z_k,   Z_k = fraction of component k inside the shell,
  %
  %   so mu, Sigma and ComponentProportion describe the real populations,
  %   including the part that cannot be seen, and temperatures and densities
  %   built from them are not biased by the cut. The properties of
  %   gmdistribution that GMDist uses are provided (mu, Sigma,
  %   ComponentProportion, NumComponents, BIC, pdf, posterior, ...).
  %
  %   Fit:  g = GMTrunc.fit(X,K,[vmin vmax],option,value,...)   X: [N x 3] km/s
  %   Test: GMTrunc.demo
  %
  %   See also GMDIST, FITGMDIST, GMDISTRIBUTION

  properties (SetAccess = protected)
    mu                    % [K x 3] means of the untruncated components [km/s]
    Sigma                 % [3 x 3 x K] covariances [(km/s)^2]
    ComponentProportion   % [1 x K] proportions of the untruncated mixture
    vRange                % [vmin vmax] [km/s]
    Z                     % [1 x K] fraction of each component inside the shell
    NumComponents
    NumVariables = 3
    CovarianceType = 'full'
    SharedCovariance = false
    nMC = 1e5             % random draws used for Z and the missing moments
    BIC = NaN
    NegativeLogLikelihood = NaN
    NumIterations = NaN
    Converged = false
  end

  methods
    function obj = GMTrunc(mu,Sigma,p,vRange,varargin)
      % GMTRUNC Construct from the parameters of the untruncated components.
      %
      %   g = GMTrunc(mu,Sigma,p,[vmin vmax][,'nMC',n])
      q = inputParser;
      addParameter(q,'nMC',1e5);
      parse(q,varargin{:});
      obj.mu = mu;
      obj.Sigma = Sigma;
      obj.ComponentProportion = p(:).'/sum(p);
      obj.vRange = vRange(:).';
      obj.NumComponents = size(mu,1);
      obj.nMC = q.Results.nMC;
      obj.Z = GMTrunc.insideFraction_(mu,Sigma,obj.vRange,obj.nMC);
    end

    function A = acceptance(obj)
      % ACCEPTANCE Fraction of the whole mixture that falls inside the shell.
      A = sum(obj.ComponentProportion.*obj.Z);
    end

    function p = pdf(obj,X)
      % PDF Density of the observed (truncated) mixture at the rows of X.
      sp = sqrt(sum(X.^2,2));
      inside = sp > obj.vRange(1) & sp < obj.vRange(2);
      p = zeros(size(X,1),1);
      d = zeros(nnz(inside),1);
      for k = 1:obj.NumComponents
        d = d + obj.ComponentProportion(k)*exp(GMTrunc.logGpdf_(X(inside,:),obj.mu(k,:),obj.Sigma(:,:,k)));
      end
      p(inside) = d/obj.acceptance();
    end

    function [P,nlogl] = posterior(obj,X)
      % POSTERIOR Responsibilities [N x K] of the components for the rows of X.
      lp = zeros(size(X,1),obj.NumComponents);
      for k = 1:obj.NumComponents
        lp(:,k) = log(obj.ComponentProportion(k)) + GMTrunc.logGpdf_(X,obj.mu(k,:),obj.Sigma(:,:,k));
      end
      mx = max(lp,[],2);
      w = exp(lp - mx);
      tot = sum(w,2);
      P = w./tot;
      nlogl = -sum(mx + log(tot)) + size(X,1)*log(obj.acceptance());
    end

    function X = random(obj,N)
      % RANDOM Draw N points from the observed (truncated) mixture.
      X = zeros(0,3);
      cs = cumsum(obj.ComponentProportion);
      for it = 1:1000
        M = max(ceil(2*N/max(obj.acceptance(),1e-3)),1000);
        idx = 1 + sum(rand(M,1) >= cs(1:end-1),2);
        Y = zeros(M,3);
        for k = 1:obj.NumComponents
          ii = find(idx == k);
          if isempty(ii), continue, end
          L = chol(obj.Sigma(:,:,k),'lower');
          Y(ii,:) = obj.mu(k,:) + randn(numel(ii),3)*L.';
        end
        sp = sqrt(sum(Y.^2,2));
        X = [X; Y(sp > obj.vRange(1) & sp < obj.vRange(2),:)]; %#ok<AGROW>
        if size(X,1) >= N, X = X(1:N,:); return, end
      end
      error('GMTrunc:random','The shell holds too little of the mixture to draw from.')
    end

    function g2 = similar(obj,mu,Sigma,p)
      % SIMILAR A mixture with other parameters but the same shell.
      g2 = GMTrunc(mu,Sigma,p,obj.vRange,'nMC',obj.nMC);
    end
  end

  methods (Static)
    function obj = fit(X,K,vRange,varargin)
      % FIT Maximum-likelihood fit of a K-component Gaussian mixture, truncated to a shell.
      %
      %   g = GMTrunc.fit(X,K,[vmin vmax][,option,value,...])
      %
      %   X: observed velocities [N x 3] (km/s), all inside the shell. EM with the
      %   missing part of each component (outside the shell) added to the sums of
      %   the E and M steps; its mass and moments come from a fixed set of random
      %   draws, so the iteration is deterministic. The log-likelihood is that of
      %   the observed points, sum log sum_k pi_k N_k - N log A, and the BIC uses it.
      %   Options:
      %     'Start'               struct with mu, Sigma, ComponentProportion (for
      %                           example the previous fit), or anything else: start
      %                           from a plain fitgmdist fit of X (default)
      %     'RegularizationValue' added to the diagonal of each Sigma [(km/s)^2] (0.1)
      %     'MaxIter'             150
      %     'Tol'                 stop when no mean moves more than this [km/s] and no
      %                           proportion more than 1e-4 (0.05)
      %     'nMC'                 random draws for the missing part (1e5)
      %   Other options (Options, SharedCovariance, ...) go to the starting fitgmdist.
      q = inputParser;
      q.KeepUnmatched = true;
      addParameter(q,'Start','gmm');
      addParameter(q,'RegularizationValue',0.1);
      addParameter(q,'MaxIter',150);
      addParameter(q,'Tol',0.05);
      addParameter(q,'nMC',1e5);
      parse(q,varargin{:});
      reg = q.Results.RegularizationValue;
      nMC = q.Results.nMC;
      vmin = vRange(1); vmax = vRange(2);
      N = size(X,1);

      st = q.Results.Start;
      if isstruct(st)
        pi_ = st.ComponentProportion(:).';
        mu = st.mu;
        S = st.Sigma;
      else
        extra = [fieldnames(q.Unmatched), struct2cell(q.Unmatched)]';
        g0 = fitgmdist(X,K,'RegularizationValue',reg,extra{:});
        pi_ = g0.ComponentProportion(:).';
        mu = g0.mu;
        S = g0.Sigma;
      end

      z = GMTrunc.fixedNormals_(nMC); % randum numbers with the same (local) seed
      converged = false;
      for iter = 1:q.Results.MaxIter
        % E step on the observed points (expectation)
        lp = zeros(N,K);
        for k = 1:K % sum(log( pi_k * N(x|mu,sigma) )) = sum log(pi)+log(N(x|mu,sigma))
          lp(:,k) = log(pi_(k)) + GMTrunc.logGpdf_(X,mu(k,:),S(:,:,k));
        end
        mx = max(lp,[],2); % max likelihood (?) component
        w = exp(lp - mx); % for normalization ?, also remove the log part, so it's only the p
        tot = sum(w,2); % sum over k
        g = w./tot; % normalize
        ll = sum(mx + log(tot)); % log likelihood?
        S0 = sum(g,1);                          % [1 x K]
        S1 = g.'*X;                             % [K x 3]
        S2 = zeros(3,3,K);
        for k = 1:K, S2(:,:,k) = (X.*g(:,k)).'*X; end

        % The missing part: what each component has outside the shell
        Zk = zeros(1,K); m0 = zeros(1,K); m1 = zeros(K,3); m2 = zeros(3,3,K);
        for k = 1:K
          L = chol(S(:,:,k),'lower');
          xs = mu(k,:) + z*L.';
          sp = sqrt(sum(xs.^2,2));
          out = ~(sp > vmin & sp < vmax); % all the MC (z) points that are not inside the given v range
          xo = xs(out,:); % total speed of particles
          % the out-of-range moments
          m0(k) = nnz(out)/nMC;   % nnz, number of non-zero elements
                                  % m0, fraction of density outside v range,
                                  % assuming the microparticles carry the
                                  % same weight
          Zk(k) = 1 - m0(k);      % fraction inside v range
          m1(k,:) = sum(xo,1)/nMC;    % bulk/average speed
          m2(:,:,k) = (xo.'*xo)/nMC;  % stress 
        end
        A = sum(pi_.*Zk);         % total density fraction inside v shell
        Nt = N/A;                 % observed + missing points, N is observed points
        S0 = S0 + Nt*pi_.*m0;     % total density?
        S1 = S1 + Nt*pi_(:).*m1;
        for k = 1:K, S2(:,:,k) = S2(:,:,k) + Nt*pi_(k)*m2(:,:,k); end
        ll = ll - N*log(A); % log likelihood of observed minus (normalized to) what is inside the shell

        % M step (maximization)
        S0 = max(S0,realmin);
        piN = S0/sum(S0);
        muN = S1./S0(:);
        SN = zeros(3,3,K);
        for k = 1:K
          Sk = S2(:,:,k)/S0(k) - muN(k,:).'*muN(k,:) + reg*eye(3);
          SN(:,:,k) = (Sk + Sk.')/2;
        end
        dmu = max(abs(muN(:) - mu(:)));
        dpi = max(abs(piN(:) - pi_(:)));
        pi_ = piN; mu = muN; S = SN;
        if dmu < q.Results.Tol && dpi < 1e-4
          converged = true;
          break
        end
      end

      obj = GMTrunc(mu,S,pi_,vRange,'nMC',nMC);
      nPar = K*3 + K*6 + (K-1);
      obj.NegativeLogLikelihood = -ll;
      obj.BIC = -2*ll + nPar*log(N);
      obj.NumIterations = iter;
      obj.Converged = converged;
    end

    function demo()
      % DEMO Synthetic check: a cold core and a hot beam cut by the upper speed limit.
      %   A plain mixture fitted to the truncated sample is biased; GMTrunc recovers
      %   the truth (about 1 minute).
      rng(3)
      vRange = [50 700];
      pT = [0.6 0.4]; muT = [-350 0 0; -650 50 0]; sT = [60 120];
      N0 = 400000;
      comp = 1 + (rand(N0,1) > pT(1));
      Xall = zeros(N0,3);
      for k = 1:2
        ii = find(comp == k);
        Xall(ii,:) = muT(k,:) + sT(k)*randn(numel(ii),3);
      end
      sp = sqrt(sum(Xall.^2,2));
      X = Xall(sp > vRange(1) & sp < vRange(2),:);
      fprintf('observed %d of %d points\n',size(X,1),N0)
      fprintf('truth:   pi = 0.60 / 0.40, mu_x = -350 / -650, sigma = 60 / 120\n')
      st = struct('ComponentProportion',[0.5 0.5],'mu',[-300 0 0; -500 0 0], ...
        'Sigma',cat(3,80^2*eye(3),80^2*eye(3)));
      gp = fitgmdist(X,2,'Start',st,'RegularizationValue',0.1);
      GMTrunc.report_('plain mixture ',gp.ComponentProportion,gp.mu,gp.Sigma)
      gt = GMTrunc.fit(X,2,vRange,'Start',st);
      GMTrunc.report_('GMTrunc       ',gt.ComponentProportion,gt.mu,gt.Sigma)
      fprintf('GMTrunc: %d iterations, converged %d, fraction inside the shell per component %s\n', ...
        gt.NumIterations,gt.Converged,mat2str(gt.Z,3))
    end
  end

  methods (Static, Access = private)
    function Z = insideFraction_(mu,Sigma,vRange,nMC)
      % Fraction of each Gaussian inside the shell, from a fixed random set.
      z = GMTrunc.fixedNormals_(nMC);
      K = size(mu,1);
      Z = zeros(1,K);
      for k = 1:K
        L = chol(Sigma(:,:,k),'lower');
        xs = mu(k,:) + z*L.';
        sp = sqrt(sum(xs.^2,2));
        Z(k) = mean(sp > vRange(1) & sp < vRange(2));
      end
    end

    function z = fixedNormals_(n)
      % Standard normal draws that do not change from call to call, and do not
      % touch the global random state.
      s = RandStream('mt19937ar','Seed',20240101);
      z = randn(s,n,3);
    end

    function lp = logGpdf_(X,mu,S)
      % Log of the Gaussian density at the rows of X.
      R = chol(S);                              % S = R'R
      y = (X - mu)/R;                           % ||y||^2 = (x-mu)' S^-1 (x-mu)
      lp = -0.5*sum(y.^2,2) - sum(log(diag(R))) - 1.5*log(2*pi);
    end

    function report_(name,p,mu,S)
      [~,o] = sort(squeeze(S(1,1,:)+S(2,2,:)+S(3,3,:)));
      for k = o(:).'
        fprintf('%s pi = %.3f  mu_x = %6.0f  sigma = %5.0f\n',name,p(k),mu(k,1),sqrt(trace(S(:,:,k))/3))
      end
    end
  end
end
