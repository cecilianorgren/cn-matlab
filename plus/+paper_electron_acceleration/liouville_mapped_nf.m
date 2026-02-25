function varargout = liouville_mapped_nf(n0,T0,v0,psimax,vpsi)

% f0
units = irf_units;

n = n0*1e6; % m^-3
T = T0; % eV
vd = v0; % ms-1
vt = sqrt(2*units.e*T./units.me); % m/s

% str_info = {'unperturbed f:';...
%             ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
%             ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
%             ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
%             };           


%% Gaussian potential
if numel(psimax) == 1
  n_psi = 2;
  psi_vec = linspace(0,psimax,n_psi);
elseif numel(psimax) > 1
  psi_vec = psimax;
end
psimax_vt = sqrt(2*units.e*psimax./units.me); % m/s

vph = vpsi*1e3; % ms-1
v_max = 200e6;
nv = 10000;
v_vec = linspace(-v_max,v_max,nv);
dv = v_vec(2)-v_vec(1);
[PSI,V] = meshgrid(psi_vec,v_vec); PSI = permute(PSI,[2 1]); V = permute(V,[2 1]);
VPH = 0*PSI + vph;

%% symbolic expression, 120 s

U = units.me*(V-vph).^2/2 - units.e*PSI;
%PHI = repmat(tocolumn(double(phi(x_vec,lx))),1,nv);
%f0 = @(v) n(1)./pi^0.5./vt(1).*exp(-(v-vd(1)).^2./vt(1).^2);
%v0 = @(x,v) vph + sign(v-vph).*((v-vph).^2-2*units.e*phi(x)./units.me).^0.5;

%% Perform liouville mapping
[Fflat,Fflat_free,Fflat_trap] = get_f_flat(V,n,vt,vd,PSI,VPH);
F_ = Fflat_free;
if vph < 0
  F_(V<vph) = NaN;
elseif vph > 0
  F_(V>vph) = NaN;
elseif vph == 0
  F_(V<vph) = NaN;
end
n_ = nansum(F_,2)*dv;
nv_ = nansum(F_.*V,2)*dv;
v_ = nv_./n_;

if nargout == 1
  varargout{2} = n_(end);
elseif nargout == 2
  varargout{1} = n_(1);
  varargout{2} = n_(end);
elseif nargout == 4
  varargout{1} = n_(1);
  varargout{2} = n_(end);
  varargout{3} = v_(1);
  varargout{4} = v_(end);
end


%plot(v_vec*1e-6,F_(1,:)',v_vec*1e-6,F_(2,:)')

if 0 % plot
  figure(19); 
  h = setup_subplots(3,1);
  isub = 1;
  hca = h(isub); isub = isub + 1;
  imagesc(hca,psi_vec,v_vec*1e-6,(F_'))
  hb = colorbar('peer',hca);
  hca.YDir = 'normal';
  %caxis(gca,[-6,-2.5])
  %pause(0.1)  
  
  hca = h(isub); isub = isub + 1;
  plot(hca,psi_vec,n_*1e-6)
  
  hca = h(isub); isub = isub + 1;
  plot(hca,psi_vec,v_*1e-6)
  1;
end
if n_(end) == 0
  1;
end
