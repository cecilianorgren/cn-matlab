classdef GMDist
  % POPFIT Gaussian-mixture description of velocity distributions, for
  %   identifying plasma populations.
  %
  %   The basic input is a {nT x nK} cell of gmdistribution objects (one
  %   column per number of components K). Particle data is an optional
  %   companion.
  %
  %   Construction:
  %     gd = GMDist(gmmCell)                      % the basics
  %     gd = GMDist(gmmCell,'time',t,'mass',m)    % metadata
  %     gd = GMDist(gmmCell,'data',mp)            % companion data (see below)
  %     gd = GMDist.fit(mp,K,'start','previous')  % run the fits
  %
  %   Works from gmmCell alone: BIC (if the fits carry it), kneeK,
  %   populations (proportion, V, T), cold, marginal, plotBIC,
  %   plotComponents, and the GMM part of plotReduced.
  %   Needs companion data: density per component (unless 'n' is given to
  %   populations), particle overlays in plotReduced, BIC for fits that
  %   carry none.
  %
  %   Companion data is a [N x 3] matrix, a {nT x 1} cell of [N x 3]
  %   (velocities in km/s), or an object with methods samples(it) and
  %   length, optionally a property time, a property mass and a method n
  %   (a time series of total density in cm^-3). MParts is such an object,
  %   but GMDist does not depend on it.
  %
  %   Velocities are in km/s and temperatures in eV.
  %
  %   See also MPARTS, FITGMDIST, GMDISTRIBUTION

  properties (SetAccess = protected)
    time    % GenericTimeArray or numeric [nT x 1], [] if unknown
    gmm     % {nT x nK} gmdistribution ([] where a fit failed)
    K       % [1 x nK] number of components of each column
    bic     % [nT x nK]
    mass    % particle mass [kg], [] if unknown
    data    % companion data, [] if none
    options % options used by GMDist.fit
  end

  methods
    function obj = GMDist(gmmCell,varargin)
      if nargin < 1 || ~iscell(gmmCell) || ndims(gmmCell) > 2
        error('GMDist:input','GMDist requires a {nT x nK} cell of gmdistribution objects.')
      end
      p = inputParser;
      addParameter(p,'time',[]);
      addParameter(p,'mass',[]);
      addParameter(p,'data',[]);
      addParameter(p,'bic',[]);
      addParameter(p,'options',struct());
      parse(p,varargin{:});

      [nT,nK] = size(gmmCell);
      isG = cellfun(@(g) isempty(g) || isa(g,'gmdistribution'),gmmCell);
      if ~all(isG(:))
        error('GMDist:input','gmmCell must hold gmdistribution objects ([] for failed fits).')
      end

      data = p.Results.data;
      if isnumeric(data) && ~isempty(data), data = {data}; end
      if ~isempty(data) && ~iscell(data), GMDist.checkCompanion_(data); end
      t = p.Results.time;
      if ~isempty(t) && GMDist.lengthOf_(t) ~= nT
        error('GMDist:input','time must have %d elements.',nT)
      end
      if ~isempty(data) && GMDist.lengthOf_(data) ~= nT
        error('GMDist:input','data must have %d time steps.',nT)
      end

      % Number of components of each column
      K = nan(1,nK);
      for iK = 1:nK
        ks = cellfun(@GMDist.numComp_,gmmCell(:,iK));
        ks = unique(ks(~isnan(ks)));
        if numel(ks) > 1
          error('GMDist:input','Column %d mixes different numbers of components.',iK)
        elseif ~isempty(ks)
          K(iK) = ks;
        end
      end

      obj.gmm = gmmCell;
      obj.K = K;
      obj.time = t;
      obj.mass = p.Results.mass;
      obj.data = data;
      obj.options = p.Results.options;

      % BIC: given, else from the fits, else recomputed from the data
      B = p.Results.bic;
      if isempty(B)
        B = nan(nT,nK);
        for iK = 1:nK
          for it = 1:nT
            g = gmmCell{it,iK};
            if isempty(g), continue, end
            b = GMDist.ownBIC_(g);
            if isnan(b) && ~isempty(data)
              b = GMDist.bicFromSamples_(g,GMDist.getSamples_(data,it));
            end
            B(it,iK) = b;
          end
        end
        missing = isnan(B) & ~cellfun(@isempty,gmmCell);
        if any(missing(:))
          warning('GMDist:bic','BIC is missing for some fits; give ''data'' or ''bic''.')
        end
      elseif ~isequal(size(B),[nT,nK])
        error('GMDist:input','bic must be [%d x %d].',nT,nK)
      end
      obj.bic = B;
    end

    function kBest = kneeK(obj,varargin)
      % KNEEK Number of components at the knee of BIC vs K, per time.
      %
      %   kBest = gd.kneeK([,'method',m])
      %     'chord' (default)  the K where BIC lies farthest below the
      %                        straight line from the first to the last K
      %     'min'              the K with minimum BIC
      p = inputParser;
      addParameter(p,'method','chord');
      parse(p,varargin{:});
      method = validatestring(p.Results.method,{'chord','min'});

      [Ks,iSort] = sort(obj.K);
      B = obj.bic(:,iSort);
      nT = size(B,1);
      kBest = nan(nT,1);
      for it = 1:nT
        ok = ~isnan(B(it,:));
        if ~any(ok), continue, end
        k = Ks(ok); b = B(it,ok);
        if strcmp(method,'min') || numel(k) < 3
          [~,i] = min(b);
        else
          x = (k-k(1))/(k(end)-k(1));
          y = (b-min(b))/max(max(b)-min(b),eps);
          d = (y(1) + (y(end)-y(1))*x) - y;
          [~,i] = max(d);
        end
        kBest(it) = k(i);
      end
    end

    function pop = populations(obj,K,varargin)
      % POPULATIONS Per-component properties for fits with K components.
      %
      %   pop = gd.populations(K,[,'sort',s][,'n',nTot])
      %
      %   Fields (components in columns, sorted by s):
      %     prop  [nT x K]        proportion
      %     n     [nT x K]        density [cm^-3] (NaN without total density)
      %     V     [nT x K x 3]    mean velocity [km/s]
      %     T     [nT x K x 3 x 3] temperature tensor [eV]
      %     Tiso  [nT x K]        trace(T)/3 [eV]
      %   'sort': 'T' (default, coldest first), 'prop' (largest first), 'none'.
      %   'n': total density [cm^-3], [nT x 1]; default from the n method of
      %   the companion data, if it has one that works (e.g. MParts with weights).
      p = inputParser;
      addParameter(p,'sort','T');
      addParameter(p,'n',[]);
      parse(p,varargin{:});
      srt = validatestring(p.Results.sort,{'T','prop','none'});

      iK = find(obj.K == K,1);
      if isempty(iK)
        error('GMDist:K','No fits with K = %d. Available: %s.',K,mat2str(obj.K))
      end
      fac = obj.massFactor_();
      nT = size(obj.gmm,1);
      nTot = obj.totalDensity_(p.Results.n);

      prop = nan(nT,K); n = nan(nT,K); Tiso = nan(nT,K);
      V = nan(nT,K,3); T = nan(nT,K,3,3);
      for it = 1:nT
        g = obj.gmm{it,iK};
        if isempty(g), continue, end
        S = GMDist.fullSigma_(g);
        pr = g.ComponentProportion(:);
        mu = g.mu;
        Tc = fac*S;                                    % 3 x 3 x K
        ti = reshape(Tc(1,1,:)+Tc(2,2,:)+Tc(3,3,:),[K,1])/3;
        switch srt
          case 'T',    [~,o] = sort(ti,'ascend');
          case 'prop', [~,o] = sort(pr,'descend');
          otherwise,   o = (1:K).';
        end
        prop(it,:) = pr(o).';
        Tiso(it,:) = ti(o).';
        n(it,:) = nTot(it)*pr(o).';
        V(it,:,:) = reshape(mu(o,:),[1 K 3]);
        for k = 1:K
          T(it,k,:,:) = reshape(Tc(:,:,o(k)),[1 1 3 3]);
        end
      end
      pop = struct('K',K,'prop',prop,'n',n,'V',V,'T',T,'Tiso',Tiso);
    end

    function [isCold,pop] = cold(obj,K,varargin)
      % COLD Flag cold components (same layout as pop.Tiso).
      %
      %   isCold = gd.cold(K,'Tmax',Tm)     Tiso below Tm [eV]
      %   isCold = gd.cold(K,'ratio',r)     Tiso below r*max(Tiso) at that time
      %   Both given: both must hold.
      p = inputParser;
      addParameter(p,'Tmax',[]);
      addParameter(p,'ratio',[]);
      parse(p,varargin{:});
      if isempty(p.Results.Tmax) && isempty(p.Results.ratio)
        error('GMDist:cold','Give ''Tmax'' and/or ''ratio''.')
      end
      pop = obj.populations(K);
      isCold = true(size(pop.Tiso));
      if ~isempty(p.Results.Tmax)
        isCold = isCold & pop.Tiso < p.Results.Tmax;
      end
      if ~isempty(p.Results.ratio)
        isCold = isCold & pop.Tiso < p.Results.ratio*max(pop.Tiso,[],2);
      end
    end

    function gm = marginal(obj,it,K,dims)
      % MARGINAL Marginal mixture over velocity dimensions dims (subset of 1:3).
      iK = find(obj.K == K,1);
      if isempty(iK), error('GMDist:K','No fits with K = %d.',K), end
      g = obj.gmm{it,iK};
      if isempty(g), error('GMDist:fit','No fit at time index %d.',it), end
      S = GMDist.fullSigma_(g);
      gm = gmdistribution(g.mu(:,dims),S(dims,dims,:),g.ComponentProportion);
    end

    function h = plotBIC(obj,ax)
      % PLOTBIC BIC - min(BIC) over time and K, with the knee overlaid.
      if nargin < 2, ax = gca; end
      [Ks,iS] = sort(obj.K);
      B = obj.bic(:,iS);
      B = B - min(B,[],2,'omitnan');
      [x,isDate] = obj.xTime_();
      h = imagesc(ax,x,Ks,B.');
      axis(ax,'xy'); ylabel(ax,'K'); title(ax,'BIC - min(BIC)'); colorbar(ax);
      hold(ax,'on'); plot(ax,x,obj.kneeK(),'w.'); hold(ax,'off')
      if isDate, datetick(ax,'x','keeplimits'); end
    end

    function ax = plotComponents(obj,K,varargin)
      % PLOTCOMPONENTS Proportion (or density), speed and temperature of
      % the components over time. Extra arguments go to populations.
      pop = obj.populations(K,varargin{:});
      [x,isDate] = obj.xTime_();
      tiledlayout(3,1);
      ax(1) = nexttile;
      if all(isnan(pop.n(:)))
        plot(x,pop.prop,'.-'); ylabel('proportion')
      else
        plot(x,pop.n,'.-'); ylabel('n [cm^{-3}]')
      end
      legend(arrayfun(@(k) sprintf('comp %d',k),1:K,'UniformOutput',false))
      ax(2) = nexttile;
      plot(x,sqrt(sum(pop.V.^2,3)),'.-'); ylabel('|V| [km/s]')
      ax(3) = nexttile;
      plot(x,pop.Tiso,'.-'); ylabel('T [eV]')
      linkaxes(ax,'x')
      if isDate, for a = ax, datetick(a,'x','keeplimits'); end, end
    end

    function ax = plotReduced(obj,it,K,dims,varargin)
      % PLOTREDUCED Reduced (marginal) view of the fit at time index it.
      %
      %   gd.plotReduced(it,K,1)       1D: vx (1), vy (2) or vz (3)
      %   gd.plotReduced(it,K,[1 2])   2D: contours of the mixture
      %   Companion data adds the particle histogram / scatter.
      %   Options: 'ax', 'maxPoints' (scatter, default 5000).
      p = inputParser;
      addParameter(p,'ax',[]);
      addParameter(p,'maxPoints',5000);
      parse(p,varargin{:});
      ax = p.Results.ax;
      if isempty(ax), ax = gca; end
      if numel(dims) > 2, error('GMDist:dims','dims must have 1 or 2 elements.'), end

      gm = obj.marginal(it,K,dims);
      S = GMDist.fullSigma_(gm);
      hasData = ~isempty(obj.data);
      if hasData
        X = GMDist.getSamples_(obj.data,it);
        X = X(:,dims);
      end
      names = {'v_x','v_y','v_z'};

      nd = numel(dims);
      lo = zeros(1,nd); hi = zeros(1,nd);
      for i = 1:nd
        sd = sqrt(max(S(i,i,:),[],3));
        lo(i) = min(gm.mu(:,i)) - 4*sd;
        hi(i) = max(gm.mu(:,i)) + 4*sd;
        if hasData, lo(i) = min(lo(i),min(X(:,i))); hi(i) = max(hi(i),max(X(:,i))); end
      end

      hold(ax,'on')
      if nd == 1
        x = linspace(lo,hi,400).';
        if hasData
          histogram(ax,X,'Normalization','pdf','EdgeColor','none','FaceAlpha',0.4);
        end
        plot(ax,x,pdf(gm,x),'k','LineWidth',1.5)
        for k = 1:K
          gk = gmdistribution(gm.mu(k,:),S(:,:,k));
          plot(ax,x,gm.ComponentProportion(k)*pdf(gk,x))
        end
        xlabel(ax,[names{dims} ' [km/s]']); ylabel(ax,'pdf')
      else
        [xg,yg] = meshgrid(linspace(lo(1),hi(1),150),linspace(lo(2),hi(2),150));
        z = reshape(pdf(gm,[xg(:) yg(:)]),size(xg));
        if hasData
          idx = randperm(size(X,1),min(size(X,1),p.Results.maxPoints));
          scatter(ax,X(idx,1),X(idx,2),4,[0.6 0.6 0.6],'filled')
        end
        contour(ax,xg,yg,z,'LineWidth',1)
        plot(ax,gm.mu(:,1),gm.mu(:,2),'r+','MarkerSize',10)
        axis(ax,'equal')
        xlabel(ax,[names{dims(1)} ' [km/s]']); ylabel(ax,[names{dims(2)} ' [km/s]'])
      end
      title(ax,sprintf('time index %d, K = %d',it,K))
      hold(ax,'off')
    end
  end

  methods (Static)
    function obj = fit(data,K,varargin)
      % GMDist.fit Fit Gaussian mixtures to velocity samples.
      %
      %   gd = GMDist.fit(data,K,[option,value,...])
      %
      %   data: one [N x 3] matrix, a {nT x 1} cell of [N x 3] (km/s), or an
      %         object with methods samples(it) and length (e.g. MParts)
      %   K:    vector of component numbers
      %   Options:
      %     'start'  'randSample' (default) or 'previous' (start from the
      %              fit at the previous time step with the same K)
      %     'time','mass'  metadata (taken from the data object if not given)
      %   All other options are passed to FITGMDIST. A failed fit gives a
      %   warning and an empty entry.
      %
      %   Note: the fit is unweighted.
      if isempty(K)
        warning('GMDist:fit','Empty K, setting K = 1.')
        K = 1;
      end
      p = inputParser;
      p.KeepUnmatched = true;
      addParameter(p,'start','randSample');
      addParameter(p,'time',[]);
      addParameter(p,'mass',[]);
      parse(p,varargin{:});
      start = validatestring(p.Results.start,{'randSample','previous'});
      fitArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
      fitArgs = fitArgs(:)';

      if isnumeric(data), data = {data}; end
      t = []; m = [];
      if iscell(data)
        nT = numel(data);
      else
        GMDist.checkCompanion_(data);
        nT = data.length();
        if isprop(data,'time'), t = data.time; end
        if isprop(data,'mass'), m = data.mass; end
      end
      if ~isempty(p.Results.time), t = p.Results.time; end
      if ~isempty(p.Results.mass), m = p.Results.mass; end

      nK = numel(K);
      out = cell(nT,nK);
      tic
      for iK = 1:nK
        fprintf('GMDist.fit: K = %d (%d/%d)\n',K(iK),iK,nK);
        for it = 1:nT
          if mod(it,max(1,round(nT/10))) == 0, fprintf('  time %d/%d\n',it,nT); end
          X = GMDist.getSamples_(data,it);
          S = 'randSample';
          if strcmp(start,'previous') && it > 1 && ~isempty(out{it-1,iK})
            S = GMDist.startFrom_(out{it-1,iK});
          end
          try
            out{it,iK} = fitgmdist(X,K(iK),'Start',S,fitArgs{:});
          catch ME
            warning('GMDist:fit','Fit failed at time index %d, K = %d: %s',it,K(iK),ME.message)
          end
        end
        toc
      end
      obj = GMDist(out,'data',data,'time',t,'mass',m, ...
        'options',struct('start',start,'fitgmdist',{fitArgs}));
    end
  end

  methods (Access = private)
    function fac = massFactor_(obj)
      % (km/s)^2 -> eV for the particle mass
      u = irf_units;
      m = obj.mass;
      if isempty(m)
        m = u.mp;
        warning('GMDist:mass','Empty mass, using proton mass.')
      end
      fac = (m/u.eV)*1e6;
    end

    function nTot = totalDensity_(obj,nOpt)
      % Total density [cm^-3] per time: given, from the companion data, or NaN.
      nT = size(obj.gmm,1);
      nTot = nan(nT,1);
      if ~isempty(nOpt)
        if numel(nOpt) ~= nT, error('GMDist:n','n must have %d elements.',nT), end
        nTot = nOpt(:);
      elseif ~isempty(obj.data) && ~iscell(obj.data) && ismethod(obj.data,'n')
        try
          ts = obj.data.n;
          nTot = ts.data(:);
        catch
          % no usable density from the companion (e.g. no weights): NaN
        end
      end
    end

    function [x,isDate] = xTime_(obj)
      % Time axis for plots (falls back to the time index).
      nT = size(obj.gmm,1);
      x = (1:nT).'; isDate = false;
      t = obj.time;
      if isa(t,'GenericTimeArray')
        try
          x = t.toDatenum; x = x(:); isDate = true;
        catch
          % keep the index
        end
      elseif ~isempty(t)
        x = t(:);
      end
    end
  end

  methods (Static, Access = private)
    function n = lengthOf_(x)
      if iscell(x) || isnumeric(x), n = numel(x);
      else, n = x.length();
      end
    end

    function checkCompanion_(data)
      if ~(ismethod(data,'samples') && ismethod(data,'length'))
        error('GMDist:data',['Data must be a [N x 3] matrix, a cell of [N x 3], ' ...
          'or an object with methods samples(it) and length (such as MParts).'])
      end
    end

    function n = numComp_(g)
      if isempty(g), n = NaN; else, n = g.NumComponents; end
    end

    function X = getSamples_(data,it)
      if iscell(data), X = data{it};
      else, X = data.samples(it);
      end
    end

    function S = fullSigma_(g)
      % Covariances as 3-D array d x d x k whatever the CovarianceType.
      k = g.NumComponents; d = g.NumVariables;
      S = zeros(d,d,k);
      for ik = 1:k
        if g.SharedCovariance, s = g.Sigma; else, s = g.Sigma(:,:,ik); end
        if strcmp(g.CovarianceType,'diagonal'), s = diag(s); end
        S(:,:,ik) = s;
      end
    end

    function b = ownBIC_(g)
      try, b = g.BIC; catch, b = NaN; end
      if isempty(b), b = NaN; end
    end

    function b = bicFromSamples_(g,X)
      % BIC = -2 logL + nPar log(N), unweighted.
      k = g.NumComponents; d = g.NumVariables;
      logL = sum(log(pdf(g,X)));
      if strcmp(g.CovarianceType,'full'), pc = d*(d+1)/2; else, pc = d; end
      if g.SharedCovariance, nc = pc; else, nc = k*pc; end
      nPar = (k-1) + k*d + nc;
      b = -2*logL + nPar*log(size(X,1));
    end

    function S = startFrom_(gm)
      % Starting struct for FITGMDIST from an earlier fit.
      if any(gm.ComponentProportion == 0)
        S = 'randSample';
      else
        S = struct('mu',gm.mu,'Sigma',gm.Sigma,'ComponentProportion',gm.ComponentProportion);
      end
    end
  end
end
