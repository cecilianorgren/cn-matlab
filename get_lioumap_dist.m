function out = get_lioumap_dist(n0_orig,T0_orig,vtr_orig,vpsi_orig)
% out = get_lioumap_dist(n0_orig,T0_orig,vtr_orig,vpsi_orig)
%

% vtr = sqrt(psi/T0);
% vpsi = vpsi/vt0;

units = irf_units;
n0 = n0_orig*1e6;
T0 = T0_orig;

vt0 = sqrt(2*units.eV*T0/units.me); % m/s
psi = vtr_orig.data.^2*T0; % eV
vpsi = vpsi_orig.data*vt0; % m/s

vd0 = 0;

timeline = vtr_orig.time;
ntimes = timeline.length;

psi_vec = psi;
vph = vpsi; % ms-1
v_max = 200e6; % ms-1
nv = 10000;
v_vec = linspace(-v_max,v_max,nv);
dv = v_vec(2)-v_vec(1);
[PSI,V] = meshgrid(psi_vec,v_vec); PSI = permute(PSI,[2 1]); V = permute(V,[2 1]);


%% symbolic expression, 120 s

%U = units.me*(V-vph).^2/2 - units.e*PSI;
%PHI = repmat(tocolumn(double(phi(x_vec,lx))),1,nv);
%f0 = @(v) n(1)./pi^0.5./vt(1).*exp(-(v-vd(1)).^2./vt(1).^2);
%v0 = @(x,v) vph + sign(v-vph).*((v-vph).^2-2*units.e*phi(x)./units.me).^0.5;

%% Perform liouville mapping
F_all = nan(size(PSI));
for itime = 1:ntimes
  if not(isnan(PSI(itime,:)))
    vph = vpsi(itime);
    VPH = 0*PSI + vph;
    [Fflat,Fflat_free,Fflat_trap] = get_f_flat(V(itime,:),n0,vt0,vd0,PSI(itime,:),VPH(itime,:));    
    F_ = Fflat_free;
    if vph < 0
      F_(V(itime,:)<vph) = NaN;
    elseif vph > 0
      F_(V(itime,:)>vph) = NaN;
    elseif vph == 0
      F_(V(itime,:)<vph) = NaN;
    end
    F_all(itime,:) = F_;
  end
end
n_ = nansum(F_,2)*dv;
nv_ = nansum(F_.*V,2)*dv;
v_ = nv_./n_;

%% Assign output
PD = PDist(timeline,F_all,'line (reduced)',v_vec*1e-3);
PD.ancillary.v_edges = [v_vec(1)-dv v_vec+dv]*1e-3;  

out = PD;
