function varargout = f_inp(input)

units = irf_units;% Physical constants
qe = 1.602e-19;
me = 9.109e-31;
mi = 1.673e-27;
mime = 1836;
eps0 = 8.854e-12;
kB = 1.381e-23;
c = 299792458;

switch input
  case 'slow_e'    
    n = [0.005 0.005]*1e6;
    T = [700 220]; T_K = T*units.eV/kB; % use parallel temperature
    m = [1 1]*me;
    q = [-1 -1]*qe; 
    vd = [2000 17000]*1e3; % m/s
  case 'slow_i'
    1;
  case 'slow'
    n = [0.04 0.06 0.01]*1e6;
    T = [700 220 3000]; T_K = T*units.eV/kB; % use parallel temperature
    m = [1 1 1836]*me;
    q = [-1 -1 1]*qe; 
    vd = [2000 17000 0]*1e3; % m/s
    %    
    n = [0.01 0.01]*1e6;
    T = [100 4000]; T_K = T*units.eV/kB; % use parallel temperature
    m = [1 1836]*me;
    q = [-1 1]*qe; 
    vd = [20000 0]*1e3; % m/s
  case 'fast_e'
    n = [0.006 0.003]*1e6;
    T = [400 400]; T_K = T*units.eV/kB; % use parallel temperature
    m = [1 1]*me;
    q = [-1 -1]*qe; 
    vd = [-2000 30000]*1e3; % m/s
  case 'fast_e_mod'
    n = [0.006 0.001]*1e6;
    T = [1000 200]; T_K = T*units.eV/kB; % use parallel temperature
    m = [1 1]*me;
    q = [-1 -1]*qe; 
    vd = [-2000 50000]*1e3; % m/s
  case 'fast_i'
     1;
  case 'fast'
    1;
end

vt = sqrt(2*qe*T./m); % m/s
wp = sqrt(n.*q.^2./(m*eps0)); % Hz
Lin = c./wp;
Ld = vt./wp/sqrt(2);

varargout{1} = n;
varargout{2} = T;
varargout{3} = m;
varargout{4} = q;
varargout{5} = vd;
varargout{6} = vt;
varargout{7} = wp;
varargout{8} = Lin;
varargout{9} = Ld;