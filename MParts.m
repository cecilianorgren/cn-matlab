classdef MParts
  % MPARTS Class representing macroparticles generated from particle
  %   distribution functions.
  %
  %   The object always holds time and velocity-space coordinates. Per-
  %   particle attributes are optional and stored explicitly in obj.attr,
  %   one field per attribute, each a struct:
  %     attr.f.data   {nT x 1} cell; element it holds the values of all
  %                   macroparticles at time it
  %     attr.f.desc   description (char)
  %     attr.f.units  units (char)
  %   Permitted attribute names (permittedAttributes_):
  %     dn   number density carried by each macroparticle
  %     f    phase-space density
  %     dV   differential velocity-space volume
  %     dc   corresponding counts
  %   Names, shapes and fields are validated on every assignment, also for
  %   mp.attr.f.units = '...' and mp.attr.f.data = {...}. A missing desc
  %   gets a default description and a missing units is ''.
  %   The methods dn(), f(), dV() and dc() return the data ([] if absent).
  %   Nothing is derived behind your back.
  %
  %   Weights used by the moments (first match wins):
  %     1. dn
  %     2. f.*dV
  %     3. none: all particles are weighted equally. V and T are then the
  %        sample mean and covariance, and n is NaN (not defined).
  %
  %   Units: the density assumes SI weights. moments() and n() warn if
  %   the units of the weights used are empty, and error if they are not
  %   the expected ones (dn: m^-3, or f: s^3/m^6 together with dV: m^3/s^3).
  %   V and T do not depend on the weight units. Velocities are assumed
  %   to be in km/s (not tracked).
  %
  %   Generation information is grouped in obj.generation:
  %     generation.parent.iDep  parent-bin indices, {iDep1, iDep2, ...}
  %     generation.sampling       how particles were initialized in the FPI
  %                             bin (linear in v, v^2 or v^3); one per object
  %
  %   Supported distribution type (obj.type):
  %    '3D-cart'  Cartesian velocity coordinates vx, vy, and vz.
  %
  %   See also PDIST, PDIST.MACROPARTICLES

  properties (Constant, Access = protected)
    permittedAttributes_ = {'dn','f','dV','dc'};
  end

  properties (Access = protected)
    %coordinateSystem_ = '';
    t_ % GenericTimeArray
  end

  properties (Dependent = true)
    time
    %coordinateSystem
  end

  properties
    type = '' % distribution type, e.g. '3D-cart'
    species
    mass
    depend % {vx vy vz}, each a {nT x 1} cell
    attr = struct() % per-particle attributes, see class help
    generation = struct('parent',struct('iDep',{{}}),'sampling','') % how the particles were generated
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
      %     'dn','f','dv','dc'  per-particle attributes: a {nT x 1} cell,
      %                         or a struct with fields data, [desc], [units]
      %     'species'           species name
      %     'mass'              particle mass [kg]
      %     'idepparent'        cell {iDep1, iDep2, ...} of parent-bin indices
      %     'sampling'            initialization method in the bin

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
          case 'species', obj.species = value;
          case 'mass',    obj.mass = value;
          case 'idepparent'
            obj.generation.parent.iDep = cellfun(@(x) reshape(x,[nT,1]),value,'UniformOutput',false);
          case 'sampling',  obj.generation.sampling = value;
          otherwise
            k = find(strcmpi(name,MParts.permittedAttributes_),1);
            if isempty(k)
              error('MParts:unknownOption','Unknown option ''%s''.',name)
            end
            if isstruct(value)
              e = value;
              e.data = reshape(e.data,[nT,1]);
            else
              e = struct('data',{reshape(value,[nT,1])});
              %e.data = reshape(value,[nT,1]);
            end
            obj.attr.(MParts.permittedAttributes_{k}) = e; % validated in set.attr
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
            a.(fn{1}).data = a.(fn{1}).data(it);
          end
          obj.attr = a;
          for ii = 1:numel(obj.depend)
            obj.depend{ii} = obj.depend{ii}(it);
          end
          if ~isempty(obj.generation.parent.iDep)
            obj.generation.parent.iDep = cellfun(@(x) x(it),obj.generation.parent.iDep,'UniformOutput',false);
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
      % Validate attribute names, entries and shapes on every assignment,
      % including mp.attr.name.data = ... and mp.attr.name.units = ... .
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
        nm = names{ii};
        e = val.(nm);
        if ~isstruct(e) || ~isfield(e,'data')
          error('MParts:attr','attr.%s must be a struct with field ''data'' (and optionally ''desc'', ''units'').',nm)
        end
        extra = setdiff(fieldnames(e),{'data','desc','units'});
        if ~isempty(extra)
          error('MParts:attr','attr.%s has unknown field(s): %s. Allowed: data, desc, units.', ...
            nm,strjoin(extra(:).',', '))
        end
        if ~iscell(e.data) || ~isequal(size(e.data),[obj.length,1])
          error('MParts:attr','attr.%s.data must be an {nT x 1} cell (nT = %d).',nm,obj.length)
        end
        if ~isfield(e,'desc'),  e.desc  = MParts.defaultDesc_(nm); end
        if ~isfield(e,'units'), e.units = ''; end
        if ~ischar(e.desc) || ~ischar(e.units)
          error('MParts:attr','attr.%s.desc and attr.%s.units must be char.',nm,nm)
        end
        val.(nm) = orderfields(e,{'data','desc','units'});
      end
      obj.attr = val;
    end

    %function value = get.coordinateSystem(obj), value = obj.coordinateSystem_; end
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

    % Attribute data ([] if the attribute is absent). Units and
    % description are in obj.attr.<name>.units / .desc.
    function out = dn(obj), out = obj.attrData_('dn'); end
    function out = f(obj),  out = obj.attrData_('f');  end
    function out = dV(obj), out = obj.attrData_('dV'); end
    function out = dc(obj), out = obj.attrData_('dc'); end

    function varargout = moments(obj)
      % MOMENTS Density [cm^-3], velocity [km/s] and temperature tensor [eV].
      %
      %   moms = MP.moments;          % struct with fields n, V, T
      %   [n,V,T] = MP.moments;
      %
      %   Without weights (see class help) n is NaN.
      obj.checkWeightUnits_();
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
      obj.checkWeightUnits_();
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

    function X = samples(obj,it)
      % SAMPLES Velocities [vx vy vz] of the macroparticles at time index it, [N x 3] (km/s).
      X = [obj.depend{1}{it}(:), obj.depend{2}{it}(:), obj.depend{3}{it}(:)];
    end

    function w = weights(obj,it)
      % WEIGHTS Weight per macroparticle at time index it, [N x 1].
      %   [] if there are no weights (all particles equally weighted).
      if obj.hasWeights_(), w = obj.weights_(it);
      else, w = [];
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
    function out = attrData_(obj,name)
      if isfield(obj.attr,name), out = obj.attr.(name).data;
      else, out = [];
      end
    end

    function checkWeightUnits_(obj)
      % Check the units of the attributes used as weights (only the
      % density depends on them). Empty units -> warning, unexpected
      % units -> error.
      if ~obj.hasWeights_(), return, end
      a = obj.attr;
      if isfield(a,'dn'), names = {'dn'}; else, names = {'f','dV'}; end
      for ii = 1:numel(names)
        nm = names{ii};
        expected = MParts.expectedUnits_(nm);
        u = strrep(strrep(a.(nm).units,' ',''),'*','');
        if isempty(u)
          warning('MParts:units','attr.%s.units is empty; assuming %s.',nm,expected{1})
        elseif ~ismember(u,expected)
          error('MParts:units','attr.%s.units = ''%s'', but the density calculation assumes %s.', ...
            nm,a.(nm).units,strjoin(expected,' or '))
        end
      end
    end

    function tf = hasWeights_(obj)
      a = obj.attr;
      tf = isfield(a,'dn') || (isfield(a,'f') && isfield(a,'dV'));
    end

    function w = weights_(obj,it)
      % Weight per macroparticle at time index it (column vector).
      a = obj.attr;
      if isfield(a,'dn')
        w = a.dn.data{it}(:);
      elseif isfield(a,'f') && isfield(a,'dV')
        w = a.f.data{it}(:).*a.dV.data{it}(:);
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
    function d = defaultDesc_(name)
      switch name
        case 'dn', d = 'density per macroparticle';
        case 'f',  d = 'phase-space density';
        case 'dV', d = 'differential velocity-space volume';
        case 'dc', d = 'counts per macroparticle';
        otherwise, d = '';
      end
    end

    function u = expectedUnits_(name)
      % Accepted unit strings (spaces and '*' removed) for the weights the
      % density calculation assumes: dn in m^-3, or f in s^3/m^6 times dV
      % in m^3/s^3. Extend the lists to accept other spellings.
      switch name
        case 'dn', u = {'m^-3','1/m^3','m-3'};
        case 'f',  u = {'s^3/m^6','s^3m^-6','s3/m6'};
        case 'dV', u = {'m^3/s^3','m^3s^-3','m3/s3'};
        otherwise, u = {};
      end
    end
  end
end
