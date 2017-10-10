classdef Distribution < TSeries
  % Particle distributions, subclass of TSeries
  
  properties (Access=protected)
%     coordinateSystem_ = '';
%     data_
%     t_ % GenericTimeArray
%     fullDim_
%     tensorOrder_ = 0;
%     tensorBasis_ = '';
  end
  
  properties (Dependent = true)
%     data
%     time
%     coordinateSystem
  end
  
  properties (SetAccess = immutable,Dependent = true)
%     tensorOrder
%     tensorBasis
  end
  
  properties (Constant = true, Hidden = true)
%     MAX_TENSOR_ORDER = 2;
%     BASIS = {'xyz','rtp','rlp','rpz','xy','rp'};
%     BASIS_NAMES = {...
%       'Cartesian','Spherical,colatitude', 'Spherical,latitude','Cylindrical',...
%       'Cartesian 2D','Polar 2D'};
  end
  
  properties (SetAccess = protected)
%     representation
  end
  
  properties
%     name = '';
%     units = '';
%     siConversion = '';
%     userData = [];
  end
  
  methods
    function obj = Distribution(t,data,varargin)
      obj@TSeries(t,data,varargin{:});       
    end
    
    function pitchangles(obj,obj1)
      %PITCHANGLES Calculate pitchangle distribution
      % Distribution.pitchangles(pitchangles)
      %   pitchangles (default) - [0:15:180] 
    end
  end
  
end