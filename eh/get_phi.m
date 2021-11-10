function [phi,phi_progressive,ancillary] = get_phi(E,vph,tint1,tint2)

% Prepare electric field
% later we need just to multiply with vph
ffilt = 30; % Hz
Etoint = E;
Etoint = Etoint.filt(ffilt,0,[],3).tlim(tint1);
intEdt = irf_integrate(Etoint,tint1(1));
minpeakdistance = 40; % the one i had before
minpeakdistance = 40;
[PKS,LOCS,W] = findpeaks(intEdt.data,'MinPeakDistance',minpeakdistance);

% Remove phi baselevel
ts_detrend_locs = irf.ts_scalar(intEdt.time([1; LOCS; end]),intEdt.data([1; LOCS; end]));
ts_phi_baselevel = ts_detrend_locs.resample(intEdt);
intEdt_detrend = intEdt-ts_phi_baselevel;

% Phase speed and potential zero level shift
phi_shift = 0; % to keep potential > 0

% Potential from observed E
phi_timeline = intEdt_detrend.time;
phi_detrend = intEdt_detrend*vph*1e-3;
phi_detrend_shift = phi_detrend + phi_shift;
phi_detrend_shift_only_pos_vals = phi_detrend_shift;
phi_detrend_shift_only_pos_vals.data(phi_detrend_shift_only_pos_vals.data<0) = 0;

phi = phi_detrend_shift_only_pos_vals.tlim(tint2);

% Find time intervals for individual EHs
% phi_detrend_shift_only_pos_vals.data<0
ind_phi_0 = find(phi.data == 0);
diff_tmp = diff(ind_phi_0);
ind_phi_0(diff_tmp==1) = [];
nPhi = numel(ind_phi_0)-1;

stdPhi = std(phi.data);
for iPhi = 1:nPhi
  ind_tmp = (ind_phi_0(iPhi)+1):(ind_phi_0(iPhi+1)-1);
  phi_tmp = phi.data(ind_tmp);
  if max(phi_tmp) < 0.2*stdPhi
    phi.data(ind_tmp) = 0;
  end
end


% put very short intervals with low phi to zero.


phi_progressive.Etoint = Etoint;
phi_progressive.intEdt = intEdt;
phi_progressive.ts_detrend_locs = ts_detrend_locs;
phi_progressive.ts_phi_baselevel = ts_phi_baselevel;
phi_progressive.intEdt_detrend = intEdt_detrend;
phi_progressive.intEdt_detrend = intEdt_detrend;
phi_progressive.phi_detrend = phi_detrend;
phi_progressive.phi_detrend_shift = phi_detrend_shift;
phi_progressive.phi_detrend_shift_only_pos_vals = phi_detrend_shift_only_pos_vals;
phi_progressive.phi = phi;

t0 = tint2(1);
tcenter = tint2(1) + 0.5*(tint2(2)-tint2(1));
t_vec = phi.time - t0;
x_vec = (phi.time - tcenter)*vph;

ancillary.tint1 = tint1;
ancillary.tint2 = tint2;
ancillary.ffilt = ffilt;
ancillary.minpeakdistance = minpeakdistance;
ancillary.PKS = PKS;
ancillary.LOCS = LOCS;
ancillary.vph = vph;
ancillary.t0 = t0;
ancillary.tcenter = tcenter;
ancillary.t_vec = t_vec;
ancillary.x_vec = x_vec;
%ancillary.ind_phi_0 = ind_phi_0;