classdef MParts
  % MPARTS Class representing macroparticles generated from particle
  %   distribution functions.
  %
  %   The object always holds time and velocity-space coordinates. Per-
  %   particle attributes are optional and stored explicitly in obj.attr,
  %   exactly as given (nothing is derived behind your back):
  %     dn   number density carried by each macroparticle
  %     f    phase-space density
  %     dV   differential velocity-space volume
  %     dc   corresponding counts
  %   Each attribute is a {nT x 1} cell; element it holds the values of all
  %   macroparticles at time it. The names are listed in
  %   permittedAttributes_; other names are rejected. The four standard
  %   attributes can also be read and set directly, e.g. mp.dn, mp.dn = x;
  %   setting an attribute to [] removes it.
  %
  %   Weights used by the moments (first match wins):
  %     1. dn
  %     2. f.*dV
  %     3. none: all particles are weighted equally. V and T are then the
  %        sample mean and covariance, and n is NaN (not defined).
  %
  %   Supported distribution type (obj.type):
  %    '3D-cart'  Cartesian velocity coordinates vx, vy, and vz.
  %
  %   See also PDIST, PDIST.MACROPARTICLES

  properties (Constant, Access = protected)
    permittedAttributes_ = {'dn','f','dV','dc'};
  end

  properties (Access = protected)
    coordinateSystem_ = '';
    t_ % GenericTimeArray
  end

  properties (Dependent = true)
    time
    coordinateSystem
    %dn % density per macroparticle (attr.dn)
    %f  % phase-space density (attr.f)
    %dV % differential velocity-space volume (attr.dV)
    %dc % counts (attr.dc)
  end

  properties
    depend  % {vx vy vz}, each a {nT x 1} cell
    attr = struct() % per-particle attributes, each a {nT x 1} cell
    parent  % struct, indices to what energy/azim/polar bin etc the particle was generated from
    v_init  % how particles were initialized in the FPI bin (linear in v, v^2 or v^3); one per object
    type = '' % distribution type, e.g. '3D-cart'
    species
    mass
    units
  end

  properties (SetAccess = protected)
    representation
  end

  methods
    function obj = MParts(time,type,varargin)
      % MPARTS Construct a macroparticle object.
      %
      %   MP = MParts(time,type,vx,vy,vz,[option,value,...])
      %
      %   Options (case insensitive):
      %     'dn'         density per macroparticle, {nT x 1} cell
      %     'f'          phase-space density, {nT x 1} cell
      %     'dV'         differential volumes, {nT x 1} cell
      %     'dc'         counts, {nT x 1} cell
      %     'species'    species name
      %     'mass'       particle mass [kg]
      %     'units'      units label
      %     'idepparent' cell {iDep1, iDep2, ...} of parent-bin indices
      %     'v_init'     initialization method in the bin

      if nargin < 2
        error('MParts:nargin','MParts requires time and type.')
      end

      obj.t_ = time;
      obj.type = type;
      nT = obj.length;
      args = varargin;

      % Required coordinates -> depend
      switch obj.type
        case '3D-cart'
          if numel(args) < 3
            error('MParts:nargin','Type ''3D-cart'' requires vx, vy and vz.')
          end
          comps = {'vx','vy','vz'};
          for ii = 1:3
            obj.depend{ii} = reshape(args{ii},[nT,1]);
            obj.representation{ii} = comps(ii);
          end
          args(1:3) = [];
        otherwise
          error('MParts:unsupportedType','Unsupported type ''%s''.',obj.type)
      end

      % Optional name-value arguments
      if mod(numel(args),2) ~= 0
        error('MParts:nargin','Options must be name-value pairs.')
      end
      while ~isempty(args)
        name = args{1}; value = args{2}; args(1:2) = [];
        switch lower(name)
          case 'dc',      obj.attr.dc = reshape(value,[nT,1]);
          case 'dn',      obj.attr.dn = reshape(value,[nT,1]);
          case 'f',       obj.attr.f  = reshape(value,[nT,1]);
          case 'dv',      obj.attr.dV = reshape(value,[nT,1]);
          case 'species', obj.species = value;
          case 'mass',    obj.mass = value;
          case 'units',   obj.units = value;
          case 'idepparent'
            obj.parent.iDep = cellfun(@(x) reshape(x,[nT,1]),value,'UniformOutput',false);
          case 'v_init',  obj.info.v_init = value;
          otherwise
            error('MParts:unknownOption','Unknown option ''%s''.',name)
        end
      end
    end

    function varargout = subsref(obj,idx)
      %SUBSREF Handle indexing. obj(it) indexes in time.
      switch idx(1).type
        case '.'
          [varargout{1:max(nargout,1)}] = builtin('subsref',obj,idx);
        case '()'
          if numel(idx(1).subs) > 1
            error('MParts:subsref','Only indexing in time, obj(it), is supported.')
          end
          it = idx(1).subs{1};
          obj.t_ = builtin('subsref',obj.t_,idx(1));
          % Slice all per-particle attributes, then assign once so that
          % the validation in set.attr sees a consistent object.
          a = obj.attr;
          for fn = fieldnames(a).'
            a.(fn{1}) = a.(fn{1})(it);
          end
          obj.attr = a;
          for ii = 1:numel(obj.depend)
            obj.depend{ii} = obj.depend{ii}(it);
          end
          if isstruct(obj.parent) && isfield(obj.parent,'iDep')
            obj.parent.iDep = cellfun(@(x) x(it),obj.parent.iDep,'UniformOutput',false);
          end
          if numel(idx) > 1
            [varargout{1:max(nargout,1)}] = builtin('subsref',obj,idx(2:end));
          else
            varargout{1} = obj;
          end
        case '{}'
          error('MParts:subsref','Not a supported subscripted reference.')
      end
    end

    function e = end(obj,~,~)
      % END Last time index, so obj(end) works.
      e = obj.length;
    end

    function obj = set.attr(obj,val)
      % Validate attribute names and shapes on every assignment,
      % including mp.attr.name = ... .
      if ~isstruct(val)
        error('MParts:attr','attr must be a struct.')
      end
      names = fieldnames(val);
      bad = setdiff(names,MParts.permittedAttributes_);
      if ~isempty(bad)
        error('MParts:attr','Unknown attribute(s): %s. Permitted: %s.', ...
          strjoin(bad(:).',', '),strjoin(MParts.permittedAttributes_,', '))
      end
      for ii = 1:numel(names)
        x = val.(names{ii});
        if ~iscell(x) || ~isequal(size(x),[obj.length,1])
          error('MParts:attr','attr.%s must be an {nT x 1} cell (nT = %d).', ...
            names{ii},obj.length)
        end
      end
      obj.attr = val;
    end

    % Accessors for the standard attributes
    function value = get.dn(obj), value = obj.getAttr_('dn'); end
    function value = get.f(obj),  value = obj.getAttr_('f');  end
    function value = get.dV(obj), value = obj.getAttr_('dV'); end
    function value = get.dc(obj), value = obj.getAttr_('dc'); end
    function obj = set.dn(obj,val), obj = obj.setAttr_('dn',val); end
    function obj = set.f(obj,val),  obj = obj.setAttr_('f',val);  end
    %function obj = set.dV(obj,val), obj = obj.setAttr_('dV',val); end
    %function obj = set.dc(obj,val), obj = obj.setAttr_('dc',val); end

    function value = get.coordinateSystem(obj), value = obj.coordinateSystem_; end
    function value = get.time(obj), value = obj.t_; end

    function l = length(obj)
      % LENGTH Number of time steps.
      if isempty(obj.t_), l = 0;
      else, l = obj.t_.length();
      end
    end

    function out = N(obj)
      % N Number of macroparticles at each time, [nT x 1].
      out = cellfun(@numel,obj.depend{1});
    end

    function varargout = moments(obj)
      % MOMENTS Density [cm^-3], velocity [km/s] and temperature tensor [eV].
      %
      %   moms = MP.moments;          % struct with fields n, V, T
      %   [n,V,T] = MP.moments;
      %
      %   Without weights (see class help) n is NaN.
      [n,V,T] = obj.computeMoments_();
      time = obj.time;
      switch nargout
        case {0,1}
          moms.n = irf.ts_scalar(time,n);
          moms.V = irf.ts_vec_xyz(time,V);
          moms.T = irf.ts_tensor_xyz(time,T);
          varargout{1} = moms;
        case 3
          varargout{1} = irf.ts_scalar(time,n);
          varargout{2} = irf.ts_vec_xyz(time,V);
          varargout{3} = irf.ts_tensor_xyz(time,T);
        otherwise
          error('MParts:moments','Use 1 or 3 output arguments.')
      end
    end

    function out = n(obj)
      % N Density [cm^-3]. Requires weights (dn, or f and dV).
      if ~obj.hasWeights_()
        error('MParts:noWeights','Density needs dn, or f and dV.')
      end
      out = irf.ts_scalar(obj.time,obj.computeMoments_());
    end

    function out = v(obj)
      % V Bulk velocity [km/s].
      [~,V] = obj.computeMoments_();
      out = irf.ts_vec_xyz(obj.time,V);
    end

    function out = T(obj)
      % T Temperature tensor [eV].
      [~,~,Tt] = obj.computeMoments_();
      out = irf.ts_tensor_xyz(obj.time,Tt);
    end

    function out = gmm(obj,K,varargin)
      % GMM Fit Gaussian mixture models to the macroparticle velocities.
      %
      %   out = MP.gmm(K,[options])
      %
      %   K is a vector of component numbers; out is {nT x numel(K)} of
      %   gmdistribution objects.
      %
      %   Options:
      %     'start'  'randSample' (default) or 'previous', which starts each
      %              fit from the fit of the previous time step (same K).
      %   All other options are passed on to FITGMDIST.
      %
      %   Note: the fit is unweighted, so it treats all macroparticles as
      %   equally important.

      if isempty(K)
        warning('MParts:gmm','Empty K, setting K = 1.')
        K = 1;
      end

      p = inputParser;
      p.KeepUnmatched = true;
      addParameter(p,'start','randSample');
      parse(p,varargin{:});
      start = validatestring(p.Results.start,{'randSample','previous'});
      fitArgs = [fieldnames(p.Unmatched), struct2cell(p.Unmatched)]';
      fitArgs = fitArgs(:)';

      nT = obj.length;
      nK = numel(K);
      out = cell(nT,nK);

      for iK = 1:nK
        for it = 1:nT
          samples = [obj.depend{1}{it}(:), obj.depend{2}{it}(:), obj.depend{3}{it}(:)];
          S = 'randSample';
          if strcmp(start,'previous') && it > 1 && ~isempty(out{it-1,iK})
            S = MParts.gmmStartFrom_(out{it-1,iK});
          end
          out{it,iK} = fitgmdist(samples,K(iK),'Start',S,fitArgs{:});
        end
      end
    end

    function out = vx(obj), out = obj.depend{1}; end
    function out = vy(obj), out = obj.depend{2}; end
    function out = vz(obj), out = obj.depend{3}; end

    % ---- Not implemented yet ----
    function TS = j(~),            error('MParts:notImplemented','j is not implemented yet.'), TS = []; end
    function MP = sph(~),          error('MParts:notImplemented','sph is not implemented yet.'), MP = []; end
    function MP = cart(~),         error('MParts:notImplemented','cart is not implemented yet.'), MP = []; end
    function MP = cyl(~),          error('MParts:notImplemented','cyl is not implemented yet.'), MP = []; end
    function MP = coordtrans(~,varargin), error('MParts:notImplemented','coordtrans is not implemented yet.'), MP = []; end
  end

  methods (Access = private)
    function value = getAttr_(obj,name)
      if isfield(obj.attr,name), value = obj.attr.(name);
      else, value = [];
      end
    end

    function obj = setAttr_(obj,name,val)
      a = obj.attr;
      if isempty(val)
        if isfield(a,name), a = rmfield(a,name); end
      else
        a.(name) = val;
      end
      obj.attr = a; % validated in set.attr
    end

    function tf = hasWeights_(obj)
      tf = ~isempty(obj.dn) || (~isempty(obj.f) && ~isempty(obj.dV));
    end

    function w = weights_(obj,it)
      % Weight per macroparticle at time index it (column vector).
      if ~isempty(obj.dn)
        w = obj.dn{it}(:);
      elseif ~isempty(obj.f) && ~isempty(obj.dV)
        w = obj.f{it}(:).*obj.dV{it}(:);
      else
        w = ones(numel(obj.depend{1}{it}),1);
      end
    end

    function [n_out,V_out,T_out] = computeMoments_(obj)
      % Shared moment computation. T (and mass lookup) only if requested.
      nT = obj.length;
      n_out = zeros(nT,1);
      V_out = zeros(nT,3);
      weighted = obj.hasWeights_();
      doT = nargout >= 3;
      if doT
        T_out = zeros(nT,3,3);
        u = irf_units;
        mass = obj.mass;
        if isempty(mass)
          mass = u.mp;
          warning('MParts:mass','Empty mass, using proton mass.')
        end
      end

      for it = 1:nT
        w = obj.weights_(it);
        n_raw = sum(w);                  % m^-3 if weighted (assumes dn in m^-3)
        if isempty(w) || n_raw == 0, continue, end
        v = [obj.depend{1}{it}(:), obj.depend{2}{it}(:), obj.depend{3}{it}(:)]; % km/s
        if weighted
          n_out(it) = n_raw*1e-6;        % cm^-3
        else
          n_out(it) = NaN;
        end
        V_out(it,:) = (w.'*v)/n_raw;     % km/s
        if doT
          dv = v - V_out(it,:);
          % <w dv_i dv_j>/n_raw*(m/e) with v in km/s -> eV
          fac = (mass/u.eV)*1e6/n_raw;
          T_out(it,:,:) = reshape((dv.'*(w.*dv))*fac,[1 3 3]);
        end
      end
    end
  end

  methods (Static, Access = private)
    function S = gmmStartFrom_(gm)
      % Starting struct for FITGMDIST from an earlier fit.
      if any(gm.ComponentProportion == 0)
        S = 'randSample';
      else
        S = struct('mu',gm.mu,'Sigma',gm.Sigma,'ComponentProportion',gm.ComponentProportion);
      end
    end
  end
end
