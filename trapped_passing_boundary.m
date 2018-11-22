function vout = trapped_passing_boundary(v,Binfty,Bloc,phi)
% vout = trapped_passing_boundary(v,Binfty,Bloc,phi)
% input/output units:
% v: km/s
% B: nT
% phi: eV
% 
% plot(vout.vpar,vout.vperp)

v = v*1e3;
Binfty = Binfty*1e-9;
Bloc = Bloc*1e-9;
%phi = units.e*phi;

units = irf_units;

v_par_infty = 0;
v_perp_infty = v;
E = units.me*v.^2/2;
mu = units.me*v_perp_infty.^2/2/Binfty;

E_perp_loc = mu*Bloc;
E_par_loc = E - E_perp_loc + units.e*phi;
v_perp_loc = sqrt(2*E_perp_loc/units.me);
v_par_loc = sqrt(2*E_par_loc/units.me);

vout.vpar = v_par_loc*1e-3;
vout.vperp = v_perp_loc*1e-3;