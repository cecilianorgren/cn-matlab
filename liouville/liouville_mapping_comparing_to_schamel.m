% Compare mapped distributions and densities using my way and Schamels way.
% If they prove to be the same, I can use his analytical expressions for 
% the density to compare to data. As it is now I only do it numerically.
%
% - My way: I shift the velocities in U with vpsi and keep f0 stationary
% - Schamel's way: He shifts f0 with vd = v-vpsi


%% Set up expression for f0, vd defined separately below
units = irf_units;
n0 = 0.035*1e6;
T0 = 350; % eV
vt = sqrt(2*units.e*T0./units.me); % m/s
       
f0 = @(v,vd) n0./pi^0.5./vt.*exp(-(v-vd).^2./vt.^2);

%% Set up (v,t) grid
lt = 0.4; % s
t0 = 1.2; % s

t_vec = 0:0.03:4;
nt = numel(t_vec);
v_max = 100e6;
nv = 10000;
v_vec = linspace(-v_max,v_max,nv);
dv = v_vec(2)-v_vec(1);
[T,V] = meshgrid(t_vec,v_vec); T = permute(T,[2 1]); V = permute(V,[2 1]);

%% Set up expression for potential
phimax = 2400; % 340

phi = @(t,lt) phimax*0.5*(1+tanh((t-t0)/lt));
phi_vec = double(phi(t_vec,lt));
PHI = repmat(tocolumn(phi_vec),1,nv);

u =  @(v,vph,t,lt) units.me*(v-vph).^2/2 - units.e*double(phi(t,lt));
U = units.me*(V-vph).^2/2 - units.e*double(phi(T,lt));

%v0 = @(t,v) vph + sign(v-vph).*((v-vph).^2-2*units.e*phi(t)./units.me).^0.5;

%% Perform liouville mapping
% My way
vd_my = 0;
vph_my = -0*17000.0001e3;
[Fflat,Fflat_free,Fflat_trap,V0] = get_f_flat(V,n,vt,vd_my,1*PHI,PHI*0+vph_my);
F_my = Fflat_free;
F_my(V<vph_my) = 0; % make it dark blue instead of white (NaN)
n_my = nansum(F_my,2)*dv;
%[n_lb,n_sep] = paper_electron_acceleration.liouville_mapped_nf(n0*1e-6,T0,0,phimax,-vph*1e-3)
%n_(end)*1e-6

% Schamel's way
vd_sch = -0*vph_my;
vph_sch = 0;
[Fflat,Fflat_free,Fflat_trap,V0] = get_f_flat(V,n,vt,vd_sch,1*PHI,PHI*0+vph_sch);
F_sch = Fflat_free;
F_sch(V<vph_sch) = 0; % make it dark blue instead of white (NaN)
n_sch = nansum(F_sch,2)*dv;

%[n_lb,n_sep] = paper_electron_acceleration.liouville_mapped_nf(n0*1e-6,T0,0,phimax,-vph*1e-3)
%n_(end)*1e-6

% Analytical expression
% using Wolfram
% integrate exp(-(a+sqrt(x))^2)
% -> indefinite integral: -e^(-(a + sqrt(x))^2) - a sqrt(π) erf(a + sqrt(x))
n_indef = @(vph,U) -exp((-(vph + sqrt(U))^2)) - vph*sqrt(pi)*erf(vph + sqrt(U));

% Chen 2002a, eq. 6
fun_np = @(phi) exp(phi).*(1-erf(sqrt(phi)));
np = 1*n0*fun_np(double(phi(t_vec,lt))/T0);


for iphi = 1:nt
  [n_lb_tmp,n_sep_tmp] = paper_electron_acceleration.liouville_mapped_nf(n0*1e-6,T0,0,phi_vec(iphi),0);
  nfun(iphi) = n_sep_tmp;
end

%% Plot results for comparison
nlim = [0.001 0.039];
vlim = 39.9*[-1 1];
flim = [-6 -2.5];

nrows = 5;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % my distribution
  hca = h(isub); isub = isub + 1;
  pcolor(hca,T,V*1e-6,F_my)
  shading(hca,'flat')
  hca.YLim = vlim;
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'v (10^3 km/s)';
end
if 1 % schamel's distribution
  hca = h(isub); isub = isub + 1;
  pcolor(hca,T,V*1e-6,F_sch)
  shading(hca,'flat')
  hca.YLim = vlim+vd_sch*1e-6;
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'v+v_\phi (10^3 km/s)';
end
if 1 % schamel's density
  hca = h(isub); isub = isub + 1;
  plot(hca,t_vec,n_my*1e-6,t_vec,n_sch*1e-6)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'n (cm^{-3})';
  if 0 % plot theoretical values
    hold(hca,'on')    
    plot(hca,t_vec,np*1e-6)
    hold(hca,'off')
  end
end
if 1 % theoreticalm density
  hca = h(isub); isub = isub + 1;
  plot(hca,t_vec,0.5*np*1e-6)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'n (...)';
end
if 1 % theoreticalm density
  hca = h(isub); isub = isub + 1;
  plot(hca,t_vec,nfun*1e-6)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'n (...)';
end
