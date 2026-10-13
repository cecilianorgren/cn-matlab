classdef MParts
  % MPARTS Class representing macroparticles generated from particle
  %   distribution functions.
  %
  %   The object stores particle weights, velocity-space coordinates,
  %   differential velocity-space volumes, and associated metadata such
  %   as species, mass, units, and parent-bin indices.
  %
  %   Per-time quantities (data, dV, dc, depend{i}) are {nT x 1} cell
  %   arrays; element it holds the values of all macroparticles at time it.
  %
  %   Supported distribution type (obj.type):
  %    '3D-cart'  Cartesian velocity coordinates vx, vy, and vz.
  %
  %   obj.dataType says what obj.data holds: 'f' (phase-space density,
  %   default) or 'dn' (f.*dV).
  %
  %   See also PDIST, PDIST.MACROPARTICLES

  properties (Access = protected)
    coordinateSystem_ = '';
    data_ % {nT x 1} cell, f or dn per macroparticle
    t_    % GenericTimeArray
  end

  properties (Dependent = true)
    data
    time
    coordinateSystem
  end

  properties
    depend         % {vx vy vz}, each a {nT x 1} cell
    dV             % {nT x 1} cell, differential velocity-space volumes
    dc             % {nT x 1} cell, corresponding counts
    parent         % struct, indices to what energy/azim/polar bin etc the particle was generated from
    v_init
    type = ''      % distribution type, e.g. '3D-cart'
    dataType = 'f' % content of data: 'f' or 'dn'
    species
    mass
    units
  end

  properties (SetAccess = protected)
    representation
  end

  methods
    function obj = MParts(time,data,units,type,varargin)
      % MPARTS Construct a macroparticle object.
      %
      %   MP = MParts(time,data,units,type,vx,vy,vz,[option,value,...])
      %
      %   Options (case insensitive):
      %     'dv'         differential volumes, same layout as data
      %     'dc'         corresponding counts, same layout as data
      %     'datatype'   'f' (default) or 'dn'
      %     'species'    species name
      %     'mass'       particle mass [kg]
      %     'units'      overrides the units input
      %     'idepparent' cell {iDep1, iDep2, ...} of parent-bin indices
      %     'v_init'     initial velocities

      if nargin < 4
        error('MParts:nargin','MParts requires time, data, units and type.')
      end

      obj.t_ = time;
      obj.type = type;
      obj.units = units;
      nT = obj.length;
      obj.data_ = reshape(data,[nT,1]);

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
          case 'dv'
            obj.dV = reshape(value,[nT,1]);
          case 'dc'
            obj.dc = reshape(value,[nT,1]);
          case 'datatype'
            obj.dataType = validatestring(value,{'f','dn'});
          case 'species'
            obj.species = value;
          case 'mass'
            obj.mass = value;
          case 'units'
            obj.units = value;
          case 'idepparent'
            obj.parent.iDep = cellfun(@(x) reshape(x,[nT,1]),value,'UniformOutput',false);
          case 'v_init'
            obj.v_init = value;
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
          obj.data_ = obj.data_(it);
          if ~isempty(obj.dV), obj.dV = obj.dV(it); end
          if ~isempty(obj.dc), obj.dc = obj.dc(it); end
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

    function value = get.coordinateSystem(obj), value = obj.coordinateSystem_; end
    function value = get.data(obj), value = obj.data_; end
    function value = get.time(obj), value = obj.t_; end

    function l = length(obj)
      % LENGTH Number of time steps.
      if isempty(obj.t_), l = 0;
      else, l = obj.t_.length();
      end
    end

    function out = N(obj)
      % N Number of macroparticles at each time, [nT x 1].
      out = cellfun(@numel,obj.data_);
    end

    function varargout = moments(obj)
      % MOMENTS Density [cm^-3], velocity [km/s] and temperature tensor [eV].
      %
      %   moms = MP.moments;          % struct with fields n, V, T
      %   [n,V,T] = MP.moments;
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
      % N Density [cm^-3].
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

    function MP = dn(obj)
      % DN Return object with data = f.*dV.
      MP = obj;
      switch obj.dataType
        case 'dn'
          return
        case 'f'
          if isempty(obj.dV), error('MParts:dn','dV is empty.'), end
          MP.data_ = cellfun(@(f,dv) f(:).*dv(:),obj.data_,obj.dV,'UniformOutput',false);
          MP.dataType = 'dn';
          MP.units = '';
      end
    end

    % ---- Not implemented yet ----
    function TS = j(~),            error('MParts:notImplemented','j is not implemented yet.'), TS = []; end
    function MP = sph(~),          error('MParts:notImplemented','sph is not implemented yet.'), MP = []; end
    function MP = cart(~),         error('MParts:notImplemented','cart is not implemented yet.'), MP = []; end
    function MP = cyl(~),          error('MParts:notImplemented','cyl is not implemented yet.'), MP = []; end
    function MP = coordtrans(~,varargin), error('MParts:notImplemented','coordtrans is not implemented yet.'), MP = []; end
  end

  methods (Access = private)
    function w = weights_(obj,it)
      % Phase-space weight per macroparticle at time index it.
      switch obj.dataType
        case 'f'
          w = obj.data_{it}(:).*obj.dV{it}(:);
        case 'dn'
          w = obj.data_{it}(:);
      end
    end

    function [n_out,V_out,T_out] = computeMoments_(obj)
      % Shared moment computation. T (and mass lookup) only if requested.
      if strcmp(obj.dataType,'f') && isempty(obj.dV)
        error('MParts:moments','dV is required to compute moments from f.')
      end
      nT = obj.length;
      n_out = zeros(nT,1);
      V_out = zeros(nT,3);
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
        n_raw = sum(w);                  % m^-3 (assumes f [s^3/m^6], dV [m^3/s^3])
        if isempty(w) || n_raw == 0, continue, end
        v = [obj.depend{1}{it}(:), obj.depend{2}{it}(:), obj.depend{3}{it}(:)]; % km/s
        n_out(it) = n_raw*1e-6;          % cm^-3
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
