function varargout = ts_get_v0(phi,vphi,v_start)
% based on get_f_flatm

units = irf_units;
if nargin > 6 && species == 1
  m = units.mp;
  q = units.e;
else
  m = units.me;
  q = -units.e;
end

% Collect input
if isa(phi,'TSeries')
  phi_data = phi.data;
  returnTs = 1;
else
  phi_data = phi;
  returnTs = 0;
end

vph = vphi*1e3; % m/s
v_start = v_start*1e3; % m/s

nt = numel(phi_data);
nv = numel(v_start);
v0 = zeros(nt,nv);

for iv = 1:nv
  if v_start(iv) > vph  
    v0(:,iv) = vph + ((v_start(iv)-vph)^2 - 2*q*phi_data/m).^0.5;
  elseif v_start(iv) < vph  
    v0(:,iv) = vph - ((v_start(iv)-vph)^2 - 2*q*phi_data/m).^0.5;  
  end
end

if returnTs
  TS = irf.ts_scalar(phi.time,v0*1e-3);
  TS.units = 'km/s';
  varargout{1} = TS;
else
  varargout{1} = v0*1e-3;
end
