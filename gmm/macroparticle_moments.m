function moms = macroparticle_moments(MP, mass, time)
% MACROPARTICLE_MOMENTS Density, bulk velocity, and temperature from macroparticles.
%
%   moms = macroparticle_moments(MP)
%   moms = macroparticle_moments(MP, mass)
%   moms = macroparticle_moments(MP, mass, time)
%
% Input
%   MP    - struct array from PDist.macroparticles, one element per time
%           (fields: vx, vy, vz, df, dv). Also accepts a single struct.
%   mass  - particle mass [kg]. Default: units.mp
%   time  - EpochTT timeline. Default: EpochTT(1:numel(MP)) placeholders
%           are not created; if omitted and numel(MP)>1, outputs are numeric.
%           Prefer passing pdist.time / times from the calling script.
%
% Output
%   moms.n  - density [cm^-3]          (scalar or TSeries if time given)
%   moms.V  - bulk velocity [km/s]     (1x3 or TSeries vec_xyz)
%   moms.T  - temperature tensor [eV]  (3x3 or TSeries tensor_xyz)
%
% Notes
%   Corrects the built-in MP.mom fields, which store flux (n*V) as "vx"
%   and pressure-like quantities as "Txx" (not divided by n).

units = irf_units;
if nargin < 2 || isempty(mass)
  mass = units.mp;
end
if nargin < 3
  time = [];
end

nT = numel(MP);
n_out = zeros(nT,1);
V_out = zeros(nT,3);
T_out = zeros(nT,3,3);

for it = 1:nT
  m = MP(it);
  w = m.df(:) .* m.dv(:);   % phase-space weight per macroparticle
  n_raw = sum(w);

  if n_raw == 0 || isempty(w)
    continue
  end

  vx = m.vx(:);
  vy = m.vy(:);
  vz = m.vz(:);

  n_out(it) = n_raw * 1e-15; % s^3/km^6 * (km/s)^3 -> cm^-3
  V_out(it,:) = [sum(w.*vx), sum(w.*vy), sum(w.*vz)] / n_raw;

  dvx = vx - V_out(it,1);
  dvy = vy - V_out(it,2);
  dvz = vz - V_out(it,3);

  % <w dv_i dv_j> / n_raw * (m/e) with v in km/s -> eV
  fac = (mass/units.eV) * 1e6 / n_raw;
  T_out(it,1,1) = sum(w.*dvx.*dvx) * fac;
  T_out(it,2,2) = sum(w.*dvy.*dvy) * fac;
  T_out(it,3,3) = sum(w.*dvz.*dvz) * fac;
  T_out(it,1,2) = sum(w.*dvx.*dvy) * fac;
  T_out(it,1,3) = sum(w.*dvx.*dvz) * fac;
  T_out(it,2,3) = sum(w.*dvy.*dvz) * fac;
  T_out(it,2,1) = T_out(it,1,2);
  T_out(it,3,1) = T_out(it,1,3);
  T_out(it,3,2) = T_out(it,2,3);
end

if ~isempty(time)
  moms.n = irf.ts_scalar(time, n_out);
  moms.V = irf.ts_vec_xyz(time, V_out);
  moms.T = irf.ts_tensor_xyz(time, T_out);
else
  moms.n = n_out;
  moms.V = V_out;
  moms.T = T_out;
end
end
