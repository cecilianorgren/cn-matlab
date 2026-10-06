classdef MParts < TSeries

  properties    
    % data is f or could also be dn = f.*dv;
    depend % {vx vy vz}    
    dv
    iDepParent % indices to what energy/azim/polar bin etc the particle was generated from    
    v_init
    type % f or dn for now
    species
    mass
  end

  methods
    function obj = MParts(v,varargin)

    end
    function N = N(obj) % number of macroparticles, could be scalar array or time series (default)
    
    end
    function MP = dn(obj)
      switch type
        case 'dn'
          MP = obj;
        case 'f'
          MP = obj;
          MP.data = MP.data.*MP.dv;
          MP.type = 'dn';
          MP.units = ''; % 
      end

    end
    function TS = n(obj)

    end
    function TS = j(obj)

    end

    function MP = sph(obj)
    end
    function MP = cart(obj)
    end
    function MP = cyl(obj)
    end

    function MP = coordtrans(obj,varargin)

    end
  end